"""Generate binary event trains for the Beggs-Plenz avalanche replication.

Three temporal structures are available, all emitting the same sparse event
format so the analysis path is identical for each:

    homogeneous  every electrode fires independently at a constant rate
    bursty       a shared telegraph drive gates otherwise independent
                 electrodes, reproducing synchronised bursting without any
                 electrode-to-electrode propagation
    critical     a true Galton-Watson branching process, used as a positive
                 control to confirm the analysis recovers real criticality

Only the critical mode contains propagation. The other two are nulls: any
avalanche structure the analysis reports for them is false by construction.

Amplitudes are the supra-threshold values, drawn from the tail of either a
truncated power law or a Gaussian. For a power law the tail above a threshold
is itself a power law with the same exponent, so the amplitude marginal keeps
its heavy tail; for a Gaussian the surviving band is narrow. That contrast is
what makes the amplitude-weighted avalanche size worth measuring.
"""

import argparse
import json
import os

import numpy as np

from simulate import sample_powerlaw

REFRACTORY_MS = 20
DEFAULT_RATE_PER_HOUR = 58_000.0

# Acklam's rational approximation to the inverse normal CDF, accurate to about
# 1e-9. Needed because the Gaussian tail above a 1-in-3700 threshold cannot be
# reached by rejection sampling in reasonable time, and scipy is not a
# dependency of this project.
_PPF_A = (
    -3.969683028665376e01,
    2.209460984245205e02,
    -2.759285104469687e02,
    1.383577518672690e02,
    -3.066479806614716e01,
    2.506628277459239e00,
)
_PPF_B = (
    -5.447609879822406e01,
    1.615858368580409e02,
    -1.556989798598866e02,
    6.680131188771972e01,
    -1.328068155288572e01,
)
_PPF_C = (
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e00,
    -2.549732539343734e00,
    4.374664141464968e00,
    2.938163982698783e00,
)
_PPF_D = (
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e00,
    3.754408661907416e00,
)
_PPF_LOW = 0.02425


def norm_ppf(p):
    """Inverse standard normal CDF, vectorised over p in (0, 1)."""
    p = np.asarray(p, dtype=float)
    out = np.empty_like(p)

    lower = p < _PPF_LOW
    upper = p > 1.0 - _PPF_LOW
    central = ~(lower | upper)

    if np.any(central):
        q = p[central] - 0.5
        r = q * q
        num = (((((_PPF_A[0] * r + _PPF_A[1]) * r + _PPF_A[2]) * r + _PPF_A[3]) * r + _PPF_A[4]) * r + _PPF_A[5]) * q
        den = ((((_PPF_B[0] * r + _PPF_B[1]) * r + _PPF_B[2]) * r + _PPF_B[3]) * r + _PPF_B[4]) * r + 1.0
        out[central] = num / den

    for mask, tail, sign in ((lower, p, 1.0), (upper, 1.0 - p, -1.0)):
        if not np.any(mask):
            continue
        q = np.sqrt(-2.0 * np.log(tail[mask]))
        num = ((((_PPF_C[0] * q + _PPF_C[1]) * q + _PPF_C[2]) * q + _PPF_C[3]) * q + _PPF_C[4]) * q + _PPF_C[5]
        den = (((_PPF_D[0] * q + _PPF_D[1]) * q + _PPF_D[2]) * q + _PPF_D[3]) * q + 1.0
        out[mask] = sign * num / den

    return out


def powerlaw_icdf(u, a, xmin, xmax):
    """Quantile function of the truncated PDF proportional to x**a."""
    u = np.asarray(u, dtype=float)
    if np.isclose(a, -1.0):
        return xmin * (xmax / xmin) ** u
    if np.isclose(a, -2.0):
        return 1.0 / (1.0 / xmin - u * (1.0 / xmin - 1.0 / xmax))
    exponent = a + 1.0
    return ((xmax**exponent - xmin**exponent) * u + xmin**exponent) ** (1.0 / exponent)


def electrode_positions(count):
    """8x8 grid with the corners removed, matching the 60 channel array."""
    side = 8
    corners = {(0, 0), (0, side - 1), (side - 1, 0), (side - 1, side - 1)}
    grid = [(r, c) for r in range(side) for c in range(side) if (r, c) not in corners]
    if count <= len(grid):
        return np.array(grid[:count], dtype=np.int16)
    extra = [(r, c) for r in range(side, side + count) for c in range(1)]
    return np.array((grid + extra)[:count], dtype=np.int16)


def amplitude_threshold(args, per_slot_probability):
    """Value a sample must exceed to register as an event."""
    quantile = 1.0 - per_slot_probability
    if args.amp_dist == "powerlaw":
        return float(powerlaw_icdf(quantile, args.a, args.xmin, args.xmax))
    return float(args.mu + args.sigma * norm_ppf(quantile))


def sample_amplitudes(rng, size, args, threshold, tail_probability):
    """Draw supra-threshold amplitudes from the tail of the chosen marginal."""
    if size == 0:
        return np.empty(0)
    if args.amp_dist == "powerlaw":
        # The tail of a truncated power law is another truncated power law
        # with the same exponent, so this is exact rather than approximate.
        return sample_powerlaw(rng, size, args.a, threshold, args.xmax)
    u = rng.uniform(1.0 - tail_probability, 1.0, size)
    return args.mu + args.sigma * norm_ppf(u)


def apply_refractory(electrodes, times, refractory_ms):
    """Greedily drop events falling within the refractory window of the last kept."""
    order = np.lexsort((times, electrodes))
    electrodes, times = electrodes[order], times[order]

    same_electrode = np.zeros(times.size, dtype=bool)
    same_electrode[1:] = electrodes[1:] == electrodes[:-1]
    gaps = np.full(times.size, np.inf)
    gaps[1:] = times[1:] - times[:-1]

    keep = np.ones(times.size, dtype=bool)
    if np.any(same_electrode & (gaps < refractory_ms)):
        previous_electrode = -1
        previous_time = 0
        for index in range(times.size):
            electrode = electrodes[index]
            time = times[index]
            if electrode != previous_electrode:
                previous_electrode = electrode
                previous_time = time
            elif time - previous_time < refractory_ms:
                keep[index] = False
            else:
                previous_time = time

    return order, keep


def generate_homogeneous(rng, args, total_ms, per_slot_probability):
    """Independent Bernoulli firing at every electrode, every millisecond."""
    expected = per_slot_probability * total_ms
    electrodes, times = [], []
    for electrode in range(args.electrodes):
        count = rng.poisson(expected)
        electrodes.append(np.full(count, electrode, dtype=np.int16))
        times.append(rng.integers(0, total_ms, count, dtype=np.int64))
    return np.concatenate(electrodes), np.concatenate(times)


def generate_bursty(rng, args, total_ms, per_slot_probability):
    """A shared telegraph drive gates conditionally independent electrodes.

    Two independent knobs add heterogeneity, both lognormal and both off by
    default. Lognormal is deliberately a broad but non-power-law spread, so any
    power law that emerges is not one we injected.

    --burst-duration-sigma spreads the burst lengths. This turns out to matter
    little, because an avalanche stops at the first empty bin and so its length
    is governed by bin occupancy rather than by how long the burst runs.

    --burst-intensity-sigma spreads the within-burst rate, which is the knob
    that drives the mixture-of-geometrics mechanism. A run inside a burst of
    occupancy lambda is geometric with scale 1/(1 - lambda) below one and
    effectively unbounded above one, so spreading lambda across that boundary
    mixes wildly different scales and is the way an independent process could
    counterfeit a power law.
    """
    array_rate_per_ms = per_slot_probability * args.electrodes
    burst_rate_per_ms = array_rate_per_ms / args.burst_duty

    period_ms = args.burst_ms / args.burst_duty
    burst_count = max(1, int(total_ms / period_ms))

    spread = args.burst_duration_sigma
    if spread > 0:
        location = np.log(args.burst_ms) - 0.5 * spread**2
        durations = rng.lognormal(location, spread, burst_count)
        durations = np.clip(np.round(durations), 1, total_ms).astype(np.int64)
    else:
        durations = np.full(burst_count, args.burst_ms, dtype=np.int64)

    intensity_spread = args.burst_intensity_sigma
    if intensity_spread > 0:
        intensity = rng.lognormal(-0.5 * intensity_spread**2, intensity_spread, burst_count)
    else:
        intensity = np.ones(burst_count)

    starts = np.sort(rng.integers(0, max(1, total_ms - durations.max()), burst_count))
    per_burst = rng.poisson(burst_rate_per_ms * intensity * durations)
    total = int(per_burst.sum())

    offsets = (rng.random(total) * np.repeat(durations, per_burst)).astype(np.int64)
    times = np.repeat(starts, per_burst) + offsets
    electrodes = rng.integers(0, args.electrodes, total).astype(np.int16)
    return electrodes, times.astype(np.int64)


def generate_critical(rng, args, total_ms, per_slot_probability):
    """Galton-Watson cascades placed on the timeline, the positive control.

    An electrode already recruited by a cascade is barred from it for the rest
    of that cascade. That follows from the 20 ms refractory period: avalanches
    are far shorter than 20 ms, so an electrode physically cannot fire twice
    within one. It also gives the process the cutoff the paper measures, since
    total size can then never exceed the electrode count.

    The bar matters for a second reason. At sigma exactly 1 an unrestricted
    critical process wanders for an unbounded time before dying, and on a small
    array with only a few bins of memory it settles into a self-sustaining
    state that produces single cascades of hundreds of thousands of events.
    Exhausting the array is what actually terminates it.
    """
    bin_ms = args.critical_bin_ms
    target_events = int(per_slot_probability * args.electrodes * total_ms)
    all_electrodes = np.arange(args.electrodes)

    cascades = []
    produced = 0
    while produced < target_events:
        active = np.array([rng.integers(args.electrodes)])
        offsets = [np.zeros(1, dtype=np.int64)]
        members = [active.copy()]
        used = active.copy()
        generation = 0

        while active.size and generation < args.max_generations:
            wanted = int(rng.poisson(args.branching_sigma, active.size).sum())
            if wanted == 0:
                break
            available = np.setdiff1d(all_electrodes, used)
            if available.size == 0:
                break
            chosen = rng.choice(available, size=min(wanted, available.size), replace=False)
            generation += 1
            offsets.append(np.full(chosen.size, generation, dtype=np.int64))
            members.append(chosen)
            used = np.concatenate((used, chosen))
            active = chosen

        cascade_offsets = np.concatenate(offsets)
        cascade_members = np.concatenate(members).astype(np.int16)
        cascades.append((cascade_offsets, cascade_members))
        produced += cascade_offsets.size

    durations = np.array([offsets.max() + 1 for offsets, _ in cascades])
    total_bins = total_ms // bin_ms
    spare = max(len(cascades), int(total_bins - durations.sum()))
    gaps = 1 + rng.poisson(max(0.0, spare / len(cascades) - 1.0), len(cascades))
    starts = np.cumsum(gaps + np.concatenate(([0], durations[:-1])))

    times, electrodes = [], []
    for start, (offsets, members) in zip(starts, cascades):
        bins = start + offsets
        times.append(bins * bin_ms + rng.integers(0, bin_ms, offsets.size))
        electrodes.append(members)

    return np.concatenate(electrodes), np.concatenate(times).astype(np.int64)


GENERATORS = {
    "homogeneous": generate_homogeneous,
    "bursty": generate_bursty,
    "critical": generate_critical,
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=sorted(GENERATORS), default="homogeneous")
    parser.add_argument("--amp-dist", choices=["powerlaw", "gaussian"], default="powerlaw")
    parser.add_argument("--hours", type=float, default=70.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="data/bp_events.npz")
    parser.add_argument("--electrodes", type=int, default=60)
    parser.add_argument("--rate-per-hour", type=float, default=DEFAULT_RATE_PER_HOUR)
    parser.add_argument("--refractory-ms", type=int, default=REFRACTORY_MS)
    parser.add_argument("--a", type=float, default=-2.0)
    parser.add_argument("--xmin", type=float, default=1.0)
    parser.add_argument("--xmax", type=float, default=1e6)
    parser.add_argument("--mu", type=float, default=0.0)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--burst-ms", type=int, default=150)
    parser.add_argument("--burst-duty", type=float, default=0.0677)
    parser.add_argument(
        "--burst-duration-sigma",
        type=float,
        default=0.0,
        help="lognormal spread of burst lengths; 0 keeps every burst the same length",
    )
    parser.add_argument(
        "--burst-intensity-sigma",
        type=float,
        default=0.0,
        help="lognormal spread of the within-burst rate; 0 keeps every burst equally intense",
    )
    parser.add_argument("--branching-sigma", type=float, default=1.0)
    parser.add_argument("--critical-bin-ms", type=int, default=4)
    parser.add_argument("--max-generations", type=int, default=10_000)
    return parser.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    total_ms = int(round(args.hours * 3_600_000))
    per_slot_probability = args.rate_per_hour / 3_600_000.0 / args.electrodes
    threshold = amplitude_threshold(args, per_slot_probability)

    electrodes, times = GENERATORS[args.mode](rng, args, total_ms, per_slot_probability)
    order, keep = apply_refractory(electrodes, times, args.refractory_ms)
    electrodes, times = electrodes[order][keep], times[order][keep]
    amplitudes = sample_amplitudes(rng, times.size, args, threshold, per_slot_probability)

    realised_ms = int(times.max()) + 1 if times.size else total_ms
    meta = {
        "mode": args.mode,
        "amp_dist": args.amp_dist,
        "electrodes": args.electrodes,
        "refractory_ms": args.refractory_ms,
        "threshold": threshold,
        "target_rate_per_hour": args.rate_per_hour,
        "requested_ms": total_ms,
        "realised_ms": realised_ms,
        "seed": args.seed,
        "critical_bin_ms": args.critical_bin_ms if args.mode == "critical" else None,
        "branching_sigma": args.branching_sigma if args.mode == "critical" else None,
        "burst_duration_sigma": args.burst_duration_sigma if args.mode == "bursty" else None,
        "burst_intensity_sigma": args.burst_intensity_sigma if args.mode == "bursty" else None,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez_compressed(
        args.out,
        electrode=electrodes,
        time_ms=times,
        amplitude=amplitudes,
        positions=electrode_positions(args.electrodes),
        meta_json=np.array(json.dumps(meta)),
    )

    dropped = keep.size - keep.sum()
    realised_rate = times.size / (realised_ms / 3_600_000.0)
    print(f"mode={args.mode} amp_dist={args.amp_dist} hours={realised_ms / 3_600_000:.2f}")
    print(f"threshold={threshold:.4f}  events={times.size}  dropped_by_refractory={dropped}")
    print(f"target rate={args.rate_per_hour:.0f}/hr  realised={realised_rate:.0f}/hr")
    if amplitudes.size:
        print(
            f"amplitude range: min={amplitudes.min():.4f} "
            f"median={np.median(amplitudes):.4f} max={amplitudes.max():.4f}"
        )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
