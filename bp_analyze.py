"""Replicate the Beggs-Plenz 2003 avalanche analysis and apply it to null data.

Their avalanche definition is purely temporal. Events are binned at width dt,
a frame is the set of electrodes active in one bin, and an avalanche is a
maximal run of consecutively occupied frames bracketed by empty ones. Spatial
adjacency plays no part: the paper reports a contiguity index of only 39
percent, meaning activity usually skips its nearest neighbours.

The free parameter is dt, and the paper sets it to IEI_avg, the mean interval
between successive events anywhere on the array. Binning at that width places
roughly one event in each bin, which is exactly the condition under which the
branching parameter of an independent process equals one. Sweeping dt is
therefore the point of this script: it separates results that reflect the data
from results that reflect the binning choice.
"""

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_BINS_MS = "1,2,3,4,5,6,8,10,12,16,20,24,32,48,64,96,128"


def load_events(path):
    """Read the sparse event file, returned sorted by time."""
    with np.load(path, allow_pickle=False) as data:
        electrode = data["electrode"]
        time_ms = data["time_ms"]
        amplitude = data["amplitude"]
        meta = json.loads(str(data["meta_json"]))

    order = np.argsort(time_ms, kind="stable")
    return electrode[order], time_ms[order], amplitude[order], meta


def average_iei(times, tmax_ms):
    """Mean inter-event interval, conditioned on being below tmax.

    The paper averages the IEI distribution only out to T_max, the point where
    the cross-correlation reaches baseline, so this is a within-burst mean
    rather than the global mean. For their cultures the two differ by roughly
    fifteen fold.
    """
    intervals = np.diff(times)
    within = intervals[intervals <= tmax_ms]
    if within.size == 0:
        return float("nan")
    return float(within.mean())


def detect_avalanches(times, amplitudes, dt):
    """Group events into maximal runs of consecutively occupied bins.

    Times must be ascending. Every metric here relies on the events of a run
    occupying one contiguous slice, so unsorted input does not raise on its own,
    it just returns quietly wrong avalanches.
    """
    if times.size > 1 and np.any(np.diff(times) < 0):
        raise ValueError("times must be sorted ascending")

    bins = times // dt
    unique_bins, first_index, counts = np.unique(bins, return_index=True, return_counts=True)

    breaks = np.flatnonzero(np.diff(unique_bins) > 1) + 1
    starts = np.concatenate(([0], breaks))
    ends = np.concatenate((breaks, [unique_bins.size]))

    # Events of one run are contiguous in the time-sorted arrays, so a run
    # maps onto a single slice and the metrics come from prefix sums.
    next_index = np.concatenate((first_index[1:], [times.size]))
    event_start = first_index[starts]
    event_end = next_index[ends - 1]

    cumulative = np.concatenate(([0.0], np.cumsum(amplitudes)))

    n_ancestors = counts[starts]
    second_bin = np.minimum(starts + 1, unique_bins.size - 1)
    n_descendants = np.where((ends - starts) >= 2, counts[second_bin], 0)

    return pd.DataFrame(
        {
            "size": event_end - event_start,
            "duration": unique_bins[ends - 1] - unique_bins[starts] + 1,
            "amplitude_size": cumulative[event_end] - cumulative[event_start],
            "n_ancestors": n_ancestors,
            "n_descendants": n_descendants,
            "start_bin": unique_bins[starts],
        }
    )


def amplitude_in_threshold_units(amplitude_size):
    """Summed amplitude rescaled so the smallest avalanche is one unit."""
    smallest = amplitude_size.min()
    return amplitude_size / smallest if smallest > 0 else amplitude_size


def branching_parameter(n_ancestors, n_descendants, n_max):
    """Descendants per ancestor across the first two bins of each avalanche.

    The single-ancestor value is the clean case of their equation (1). The
    multi-ancestor value follows equations (2) and (3): descendants per
    ancestor are rounded, then weighted by ancestor count and scaled by their
    availability correction, which compensates for electrodes that are already
    occupied and so cannot be recruited.
    """
    single = n_ancestors == 1
    sigma_single = float(n_descendants[single].mean()) if np.any(single) else float("nan")

    usable = n_ancestors < n_max
    if not np.any(usable):
        return sigma_single, float("nan")

    ancestors = n_ancestors[usable].astype(float)
    descendants = np.round(n_descendants[usable] / ancestors)
    correction = (n_max - 1) / (n_max - ancestors)
    sigma_multi = float((ancestors * descendants * correction).sum() / ancestors.sum())

    return sigma_single, sigma_multi


def fit_regression(values, fit_min, fit_max):
    """Slope of a straight line through the log-log histogram, as the paper did."""
    selected = values[(values >= fit_min) & (values <= fit_max)]
    if selected.size < 10:
        return float("nan"), float("nan")

    sizes, counts = np.unique(selected, return_counts=True)
    if sizes.size < 3:
        return float("nan"), float("nan")

    log_size = np.log10(sizes)
    log_count = np.log10(counts / counts.sum())
    slope, intercept = np.polyfit(log_size, log_count, 1)

    predicted = slope * log_size + intercept
    residual = np.sum((log_count - predicted) ** 2)
    total = np.sum((log_count - log_count.mean()) ** 2)
    r_squared = 1.0 - residual / total if total > 0 else float("nan")

    return float(-slope), float(r_squared)


def fit_regression_logbins(values, bins_per_decade=10):
    """Log-log slope using logarithmically spaced bins.

    Summed amplitude spans several decades, where a per-value histogram leaves
    the whole upper tail in bins of count one and lets that noise dominate an
    unweighted regression. Logarithmic bins normalised by bin width give the
    density that the exponent actually describes.
    """
    values = np.asarray(values, dtype=float)
    values = values[values > 0]
    if values.size < 50:
        return float("nan"), float("nan")

    low, high = values.min(), values.max()
    if not high > low:
        return float("nan"), float("nan")

    n_bins = max(3, int(np.log10(high / low) * bins_per_decade))
    edges = np.logspace(np.log10(low), np.log10(high), n_bins + 1)
    counts, _ = np.histogram(values, bins=edges)

    occupied = counts > 0
    if occupied.sum() < 3:
        return float("nan"), float("nan")

    centres = np.sqrt(edges[:-1] * edges[1:])
    density = counts / (np.diff(edges) * counts.sum())

    log_centre = np.log10(centres[occupied])
    log_density = np.log10(density[occupied])
    slope, intercept = np.polyfit(log_centre, log_density, 1)

    predicted = slope * log_centre + intercept
    total = np.sum((log_density - log_density.mean()) ** 2)
    r_squared = 1.0 - np.sum((log_density - predicted) ** 2) / total if total > 0 else float("nan")

    return float(-slope), float(r_squared)


def _maximise(objective, low, high, iterations=200):
    """Golden-section search for the maximum of a unimodal objective."""
    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = low, high
    c, d = b - invphi * (b - a), a + invphi * (b - a)
    for _ in range(iterations):
        if objective(c) > objective(d):
            b, d = d, c
            c = b - invphi * (b - a)
        else:
            a, c = c, d
            d = a + invphi * (b - a)
        if abs(b - a) < 1e-8:
            break
    return 0.5 * (a + b)


def fit_mle(values, fit_min, fit_max):
    """Discrete MLE on a truncated support, compared against an exponential.

    The support is bounded because the array itself is bounded, so the
    normalising constant is a finite sum rather than a Hurwitz zeta. The
    Vuong statistic reports which of the two models fits better: positive
    favours the power law, and magnitudes below about 2 are inconclusive.
    """
    selected = values[(values >= fit_min) & (values <= fit_max)].astype(float)
    if selected.size < 50:
        return float("nan"), float("nan"), float("nan")

    support = np.arange(fit_min, fit_max + 1, dtype=float)

    # Both likelihoods depend on the data only through these two totals, so the
    # search costs one pass over the support per step rather than one pass over
    # every observation.
    count = selected.size
    total_log = float(np.log(selected).sum())
    total = float(selected.sum())

    def powerlaw_loglike(exponent):
        return -count * np.log(np.sum(support ** (-exponent))) - exponent * total_log

    def exponential_loglike(decay):
        return -count * np.log(np.sum(np.exp(-decay * support))) - decay * total

    # The lower bound is negative rather than 1, because on a bounded support
    # any exponent normalises and clipping at 1 would silently pin the estimate
    # to the boundary for flat distributions.
    tau = _maximise(powerlaw_loglike, -1.0, 6.0)
    rate = _maximise(exponential_loglike, 1e-6, 5.0)

    powerlaw_points = -np.log(np.sum(support ** (-tau))) - tau * np.log(selected)
    exponential_points = -np.log(np.sum(np.exp(-rate * support))) - rate * selected

    difference = powerlaw_points - exponential_points
    spread = difference.std()
    vuong = float(difference.sum() / (np.sqrt(difference.size) * spread)) if spread > 0 else float("nan")

    return float(tau), float(difference.sum()), vuong


def apply_control(control, electrodes, times, rng, jitter_ms, span_ms):
    """Surrogate data that destroys temporal structure while preserving rate."""
    if control == "none":
        return electrodes, times
    if control == "jitter":
        moved = times + rng.integers(-jitter_ms, jitter_ms + 1, times.size)
        moved = np.clip(moved, 0, span_ms - 1)
    elif control == "shuffle":
        moved = rng.integers(0, span_ms, times.size)
    else:
        raise ValueError(f"unknown control {control}")

    order = np.argsort(moved, kind="stable")
    return electrodes[order], moved[order]


def bin_sweep(times, amplitudes, widths, n_max, fit_max):
    """Every metric as a function of the one free analysis parameter.

    Sizes are fitted over [1, fit_max], which stops short of the electrode
    count on purpose. A bounded array puts a spike of probability exactly at
    its own size limit, because every cascade that would have grown past it is
    truncated to land there, and that spike drags an unrestricted fit well away
    from the true scaling exponent. The fit spanning the whole range is kept
    alongside as tau_regression_to_cutoff so the distortion stays visible.
    """
    rows = []
    for dt in widths:
        avalanches = detect_avalanches(times, amplitudes, dt)
        sigma_single, sigma_multi = branching_parameter(
            avalanches["n_ancestors"].to_numpy(),
            avalanches["n_descendants"].to_numpy(),
            n_max,
        )
        sizes = avalanches["size"].to_numpy()
        tau_regression, r_squared = fit_regression(sizes, 1, fit_max)
        tau_mle, llr, vuong = fit_mle(sizes, 1, fit_max)
        tau_full, r_squared_full = fit_regression(sizes, 1, n_max)
        tau_amplitude, _ = fit_regression_logbins(
            amplitude_in_threshold_units(avalanches["amplitude_size"].to_numpy())
        )

        rows.append(
            {
                "dt_ms": dt,
                "n_avalanches": len(avalanches),
                "sigma_single": sigma_single,
                "sigma_multi": sigma_multi,
                "tau_regression": tau_regression,
                "tau_regression_r2": r_squared,
                "tau_mle": tau_mle,
                "tau_regression_to_cutoff": tau_full,
                "tau_regression_to_cutoff_r2": r_squared_full,
                "powerlaw_vs_exponential_llr": llr,
                "vuong_z": vuong,
                "tau_amplitude": tau_amplitude,
                "mean_size": float(avalanches["size"].mean()),
                "max_size": int(avalanches["size"].max()),
                "mean_duration": float(avalanches["duration"].mean()),
            }
        )
    return pd.DataFrame(rows)


def plot_distribution(values, xlabel, title, path, reference_slope=None):
    sizes, counts = np.unique(values, return_counts=True)
    probability = counts / counts.sum()

    figure, axes = plt.subplots(figsize=(6, 4.5))
    axes.scatter(sizes, probability, s=26, color="#2b6cb0", label="observed")

    if reference_slope is not None and sizes.size:
        anchor = probability[0] * (sizes / sizes[0]) ** (-reference_slope)
        axes.plot(sizes, anchor, "--", color="#c53030", linewidth=1.2,
                  label=f"slope -{reference_slope}")
        axes.legend(frameon=False, fontsize=8)

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlabel(xlabel)
    axes.set_ylabel("probability")
    axes.set_title(title)
    axes.grid(True, which="both", linewidth=0.3, alpha=0.5)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def plot_phase(sweep, iei_avg, path):
    """The (sigma, tau) trajectory as a function of dt, replicating Figure 7D."""
    figure, axes = plt.subplots(figsize=(6, 4.5))
    axes.plot(sweep["sigma_single"], sweep["tau_regression"], "-o", color="#2b6cb0",
              markersize=4, linewidth=1)

    for _, row in sweep.iterrows():
        axes.annotate(f"{row['dt_ms']:.0f}", (row["sigma_single"], row["tau_regression"]),
                      fontsize=7, alpha=0.7, xytext=(3, 3), textcoords="offset points")

    axes.scatter([1.0], [1.5], marker="*", s=220, color="#c53030", zorder=5,
                 label="critical point (1, 1.5)")
    axes.axvline(1.0, linewidth=0.6, alpha=0.4)
    axes.axhline(1.5, linewidth=0.6, alpha=0.4)
    axes.set_xlabel("branching parameter sigma")
    axes.set_ylabel("size exponent tau")
    axes.set_title(f"Phase trajectory over bin width (IEI_avg = {iei_avg:.2f} ms)")
    axes.legend(frameon=False, fontsize=8)
    axes.grid(True, linewidth=0.3, alpha=0.5)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def plot_sweep(sweep, iei_avg, path):
    figure, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

    axes[0].plot(sweep["dt_ms"], sweep["sigma_single"], "-o", markersize=4, label="single ancestor")
    axes[0].plot(sweep["dt_ms"], sweep["sigma_multi"], "-s", markersize=4, label="multi ancestor")
    axes[0].axhline(1.0, color="#c53030", linestyle="--", linewidth=1)
    axes[0].set_ylabel("sigma")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(sweep["dt_ms"], sweep["tau_regression"], "-o", markersize=4, label="regression")
    axes[1].plot(sweep["dt_ms"], sweep["tau_mle"], "-s", markersize=4, label="MLE")
    axes[1].axhline(1.5, color="#c53030", linestyle="--", linewidth=1)
    axes[1].set_ylabel("tau")
    axes[1].set_xlabel("bin width dt (ms)")
    axes[1].legend(frameon=False, fontsize=8)

    for axis in axes:
        axis.axvline(iei_avg, color="#2f855a", linestyle=":", linewidth=1.2)
        axis.set_xscale("log")
        axis.grid(True, which="both", linewidth=0.3, alpha=0.5)

    axes[0].set_title(f"Metrics vs bin width (dotted line = IEI_avg = {iei_avg:.2f} ms)")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True)
    parser.add_argument("--plot-dir", required=True)
    parser.add_argument("--tmax-ms", type=float, default=200.0)
    parser.add_argument("--bins-ms", default=DEFAULT_BINS_MS)
    parser.add_argument("--control", choices=["none", "jitter", "shuffle"], default="none")
    parser.add_argument("--jitter-ms", type=int, default=4)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--native-bin-ms",
        type=int,
        default=0,
        help="summarise at this bin width instead of the IEI_avg derived one",
    )
    parser.add_argument(
        "--fit-max",
        type=int,
        default=0,
        help="largest size included in the exponent fit; 0 uses half the electrode count",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    electrodes, times, amplitudes, meta = load_events(args.events)
    span_ms = int(times.max()) + 1
    n_max = int(meta["electrodes"])

    if args.control != "none":
        electrodes, times = apply_control(
            args.control, electrodes, times, rng, args.jitter_ms, span_ms
        )
        amplitudes = rng.permutation(amplitudes)

    iei_avg = average_iei(times, args.tmax_ms)
    widths = sorted({int(round(float(w))) for w in args.bins_ms.split(",") if w.strip()})
    native = args.native_bin_ms if args.native_bin_ms > 0 else max(1, int(round(iei_avg)))
    if native not in widths:
        widths = sorted(widths + [native])

    fit_max = args.fit_max if args.fit_max > 0 else max(4, n_max // 2)
    sweep = bin_sweep(times, amplitudes, widths, n_max, fit_max=fit_max)
    sweep.insert(0, "iei_avg_ms", iei_avg)

    at_native = detect_avalanches(times, amplitudes, native)
    sweep_path = os.path.join(args.plot_dir, "bin_sweep.csv")
    sweep.to_csv(sweep_path, index=False)
    at_native.to_csv(os.path.join(args.plot_dir, "avalanches_at_iei.csv"), index=False)

    plot_distribution(at_native["size"].to_numpy(), "avalanche size (electrode activations)",
                      f"Size distribution at dt = {native} ms",
                      os.path.join(args.plot_dir, "size_distribution.png"), reference_slope=1.5)
    plot_distribution(at_native["duration"].to_numpy(), "avalanche duration (bins)",
                      f"Duration distribution at dt = {native} ms",
                      os.path.join(args.plot_dir, "duration_distribution.png"), reference_slope=2.0)
    scaled = np.round(amplitude_in_threshold_units(at_native["amplitude_size"].to_numpy()))
    plot_distribution(scaled, "avalanche size (summed amplitude, threshold units)",
                      f"Amplitude-weighted size at dt = {native} ms",
                      os.path.join(args.plot_dir, "amplitude_size_distribution.png"),
                      reference_slope=1.5)
    plot_phase(sweep, iei_avg, os.path.join(args.plot_dir, "phase_trajectory.png"))
    plot_sweep(sweep, iei_avg, os.path.join(args.plot_dir, "bin_sweep.png"))

    row = sweep[sweep["dt_ms"] == native].iloc[0]
    print(f"source: {args.events}  mode={meta['mode']}  control={args.control}")
    print(f"events={times.size}  electrodes={n_max}  hours={span_ms / 3_600_000:.2f}")
    print(f"IEI_avg (T_max={args.tmax_ms:.0f} ms) = {iei_avg:.3f} ms  -> native bin {native} ms")
    print(f"at native bin: avalanches={int(row['n_avalanches'])}  max_size={int(row['max_size'])}")
    print(f"  sigma single={row['sigma_single']:.4f}  multi={row['sigma_multi']:.4f}")
    print(f"  tau over sizes 1-{fit_max}: regression={row['tau_regression']:.4f}"
          f" (R2={row['tau_regression_r2']:.4f})  MLE={row['tau_mle']:.4f}")
    print(f"  tau over sizes 1-{n_max} (spans the cutoff)"
          f": regression={row['tau_regression_to_cutoff']:.4f}")
    print(f"  power law vs exponential: Vuong z={row['vuong_z']:.2f}"
          f"  (positive favours power law)")
    print(f"  tau amplitude-weighted={row['tau_amplitude']:.4f}")
    print(f"wrote {sweep_path}")


if __name__ == "__main__":
    main()
