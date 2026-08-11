"""Validation suite for the Beggs-Plenz replication pipeline.

Nothing in the results is worth reading until these pass. The checks fall into
three groups: exact behaviour on hand-built event trains where the right answer
can be counted by eye, recovery of known exponents from synthetic samples, and
two end-to-end checks on generated data, one against an analytical prediction
and one against the known ground truth of the positive control.

Run with: python bp_validate.py
"""

import argparse

import numpy as np

from bp_analyze import (
    branching_parameter,
    detect_avalanches,
    fit_mle,
    fit_regression,
)
from bp_simulate import (
    GENERATORS,
    apply_refractory,
    norm_ppf,
    powerlaw_icdf,
)

FAILURES = []


def check(name, condition, detail=""):
    status = "pass" if condition else "FAIL"
    print(f"  [{status}] {name}{'  ' + detail if detail else ''}")
    if not condition:
        FAILURES.append(name)


def close(a, b, tolerance):
    return abs(a - b) <= tolerance


def check_samplers():
    print("\nSamplers")
    check("norm_ppf median", close(float(norm_ppf(0.5)), 0.0, 1e-12))
    check("norm_ppf 97.5th", close(float(norm_ppf(0.975)), 1.959964, 1e-5))
    check("norm_ppf 0.1th", close(float(norm_ppf(0.001)), -3.090232, 1e-5))
    check("powerlaw icdf median a=-2", close(float(powerlaw_icdf(0.5, -2, 1, 100)), 1.980198, 1e-5))
    check("powerlaw icdf endpoints",
          close(float(powerlaw_icdf(0.0, -2, 1, 100)), 1.0, 1e-9)
          and close(float(powerlaw_icdf(1.0, -2, 1, 100)), 100.0, 1e-6))

    # a = -1 is the log-uniform special case handled separately in the icdf.
    check("powerlaw icdf median a=-1", close(float(powerlaw_icdf(0.5, -1, 1, 100)), 10.0, 1e-6))


def check_refractory():
    print("\nRefractory period")

    # Greedy semantics: 10 is dropped because it falls inside the window opened
    # by 0, but 25 survives because the window is measured from the last kept
    # event, not from the last event seen.
    electrodes = np.array([0, 0, 0], dtype=np.int16)
    times = np.array([0, 10, 25], dtype=np.int64)
    order, keep = apply_refractory(electrodes, times, 20)
    check("greedy keeps 0 and 25", list(times[order][keep]) == [0, 25], str(list(times[order][keep])))

    # The window is per electrode, so simultaneous events on different
    # electrodes must all survive.
    electrodes = np.array([0, 1, 2], dtype=np.int16)
    times = np.array([5, 5, 5], dtype=np.int64)
    order, keep = apply_refractory(electrodes, times, 20)
    check("window is per electrode", keep.sum() == 3, f"kept {keep.sum()} of 3")

    # Exactly at the boundary the event is kept, since the condition is a
    # strict inequality.
    electrodes = np.array([0, 0], dtype=np.int16)
    times = np.array([0, 20], dtype=np.int64)
    order, keep = apply_refractory(electrodes, times, 20)
    check("boundary at exactly 20 ms is kept", keep.sum() == 2)


def check_detection():
    print("\nAvalanche detection on hand-built trains")

    # Occupied milliseconds 0,1,2 then 5 then 7,8 at dt = 1.
    times = np.array([0, 0, 1, 2, 2, 2, 5, 7, 7, 8], dtype=np.int64)
    amplitudes = np.ones(times.size)
    found = detect_avalanches(times, amplitudes, 1)

    check("three avalanches", len(found) == 3, str(len(found)))
    check("sizes", list(found["size"]) == [6, 1, 3], str(list(found["size"])))
    check("durations", list(found["duration"]) == [3, 1, 2], str(list(found["duration"])))
    check("ancestors", list(found["n_ancestors"]) == [2, 1, 2], str(list(found["n_ancestors"])))
    check("descendants", list(found["n_descendants"]) == [1, 0, 1], str(list(found["n_descendants"])))

    # A single occupied frame bracketed by empty ones is a real avalanche; the
    # paper reports thousands of them.
    check("single frame counts as an avalanche", int((found["duration"] == 1).sum()) == 1)

    # Size counts repeat activations, which is why sizes can exceed the
    # electrode count in the paper.
    repeats = detect_avalanches(np.array([0, 0, 0], dtype=np.int64), np.ones(3), 1)
    check("size counts repeat activations", int(repeats["size"].iloc[0]) == 3)

    check("one empty bin splits", len(detect_avalanches(np.array([0, 2], dtype=np.int64), np.ones(2), 1)) == 2)
    check("adjacent bins merge", len(detect_avalanches(np.array([0, 1], dtype=np.int64), np.ones(2), 1)) == 1)
    check("wider bin merges the split", len(detect_avalanches(np.array([0, 2], dtype=np.int64), np.ones(2), 2)) == 1)

    summed = detect_avalanches(np.array([0, 0, 1], dtype=np.int64), np.array([2.0, 3.0, 5.0]), 1)
    check("amplitude_size sums the run", close(float(summed["amplitude_size"].iloc[0]), 10.0, 1e-12))

    # Every event must be accounted for exactly once across all avalanches.
    rng = np.random.default_rng(1)
    random_times = np.sort(rng.integers(0, 5_000, 4_000).astype(np.int64))
    for dt in (1, 3, 7):
        partition = detect_avalanches(random_times, np.ones(random_times.size), dt)
        check(f"events conserved at dt={dt}", int(partition["size"].sum()) == random_times.size)


def check_branching():
    print("\nBranching parameter")
    ancestors = np.array([1, 1, 2])
    descendants = np.array([2, 0, 4])
    single, multi = branching_parameter(ancestors, descendants, 60)

    check("single ancestor mean", close(single, 1.0, 1e-12), f"{single:.4f}")

    # Ancestor-weighted mean of rounded descendants-per-ancestor, each scaled
    # by the availability correction (nmax-1)/(nmax-na).
    expected = (1 * 2 * (59 / 59) + 1 * 0 * (59 / 59) + 2 * 2 * (59 / 58)) / 4
    check("multi ancestor with correction", close(multi, expected, 1e-9),
          f"{multi:.4f} vs {expected:.4f}")

    # An avalanche occupying every electrode has no room to recruit, so the
    # correction diverges and it must be excluded rather than produce infinity.
    _, saturated = branching_parameter(np.array([60, 1]), np.array([0, 1]), 60)
    check("saturated avalanche excluded", np.isfinite(saturated), f"{saturated:.4f}")


def check_fitting():
    print("\nExponent fitting")
    rng = np.random.default_rng(0)
    support = np.arange(1, 61, dtype=float)

    for true_tau in (1.5, 2.0, 2.5):
        pmf = support ** (-true_tau)
        sample = rng.choice(support, size=200_000, p=pmf / pmf.sum()).astype(int)
        tau, _, vuong = fit_mle(sample, 1, 60)
        check(f"MLE recovers tau={true_tau}", close(tau, true_tau, 0.03), f"got {tau:.4f}")
        check(f"Vuong favours power law at tau={true_tau}", vuong > 2, f"z={vuong:.1f}")

    for rate in (0.1, 0.3):
        pmf = np.exp(-rate * support)
        sample = rng.choice(support, size=200_000, p=pmf / pmf.sum()).astype(int)
        _, _, vuong = fit_mle(sample, 1, 60)
        regression, _ = fit_regression(sample, 1, 60)
        check(f"Vuong rejects power law for exponential rate={rate}", vuong < -2, f"z={vuong:.1f}")
        # The point of carrying both methods: regression still reports a
        # plausible-looking exponent for data that is definitely not a power law.
        check(f"regression still reports an exponent at rate={rate}",
              np.isfinite(regression), f"tau={regression:.3f}")


def build_args(**overrides):
    defaults = dict(
        mode="homogeneous", amp_dist="powerlaw", electrodes=60,
        rate_per_hour=58_000.0, refractory_ms=20, a=-2.0, xmin=1.0, xmax=1e6,
        mu=0.0, sigma=1.0, burst_ms=150, burst_duty=0.0677,
        burst_duration_sigma=0.0, burst_intensity_sigma=0.0,
        branching_sigma=1.0, critical_bin_ms=4, max_generations=10_000,
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def generate(mode, hours, seed, **overrides):
    args = build_args(mode=mode, **overrides)
    rng = np.random.default_rng(seed)
    total_ms = int(hours * 3_600_000)
    per_slot = args.rate_per_hour / 3_600_000.0 / args.electrodes
    electrodes, times = GENERATORS[mode](rng, args, total_ms, per_slot)
    order, keep = apply_refractory(electrodes, times, args.refractory_ms)
    electrodes, times = electrodes[order][keep], times[order][keep]

    # apply_refractory groups by electrode; the detector needs time order.
    by_time = np.argsort(times, kind="stable")
    return electrodes[by_time], times[by_time]


def check_homogeneous_prediction():
    print("\nAnalytical prediction: independent data gives sigma = rate x dt")
    _, times = generate("homogeneous", hours=6, seed=3)
    span_ms = int(times.max()) + 1
    rate_per_ms = times.size / span_ms

    for dt in (20, 40, 62, 80):
        found = detect_avalanches(times, np.ones(times.size), dt)
        single, _ = branching_parameter(
            found["n_ancestors"].to_numpy(), found["n_descendants"].to_numpy(), 60
        )
        predicted = rate_per_ms * dt
        check(f"sigma matches rate x dt at dt={dt} ms", close(single, predicted, 0.05),
              f"measured {single:.4f} vs predicted {predicted:.4f}")

    # The whole point: sigma = 1 is reachable on data with no propagation at
    # all, purely by choosing the bin width equal to the mean interval.
    critical_dt = int(round(1.0 / rate_per_ms))
    found = detect_avalanches(times, np.ones(times.size), critical_dt)
    single, _ = branching_parameter(
        found["n_ancestors"].to_numpy(), found["n_descendants"].to_numpy(), 60
    )
    check(f"independent data reaches sigma=1 at dt={critical_dt} ms",
          close(single, 1.0, 0.05), f"sigma={single:.4f}")

    # But it must not look like a power law, otherwise sigma and tau would both
    # be artifacts and the paper would have nothing left.
    _, _, vuong = fit_mle(found["size"].to_numpy(), 1, 60)
    check("independent data is still not a power law", vuong < -2, f"Vuong z={vuong:.1f}")


def check_critical_control():
    print("\nPositive control: a true branching process must be recovered")
    electrodes, times = generate("critical", hours=3, seed=5)
    amplitudes = np.ones(times.size)

    found = detect_avalanches(times, amplitudes, 4)
    single, _ = branching_parameter(
        found["n_ancestors"].to_numpy(), found["n_descendants"].to_numpy(), 60
    )
    sizes = found["size"].to_numpy()
    regression, _ = fit_regression(sizes, 1, 30)
    mle, _, vuong = fit_mle(sizes, 1, 30)

    check("critical sigma near 1", close(single, 1.0, 0.15), f"{single:.4f}")
    check("critical tau near 1.5 by regression", close(regression, 1.5, 0.1), f"{regression:.4f}")
    check("critical tau near 1.5 by MLE", close(mle, 1.5, 0.1), f"{mle:.4f}")
    check("critical data reads as a power law", vuong > 2, f"Vuong z={vuong:.1f}")

    # A bounded array leaves a spike of probability at exactly its own size
    # limit, so a fit that spans the cutoff reads materially lower than the
    # true exponent. This is why the reported fit stops short of it.
    spanning, _ = fit_regression(sizes, 1, 60)
    check("fitting across the cutoff biases tau downward", spanning < regression - 0.02,
          f"{spanning:.4f} spanning vs {regression:.4f} below cutoff")

    # Destroying the temporal structure must destroy the signature, otherwise
    # the pipeline is reporting something other than propagation.
    rng = np.random.default_rng(11)
    span_ms = int(times.max()) + 1
    shuffled = np.sort(rng.integers(0, span_ms, times.size))
    surrogate = detect_avalanches(shuffled, np.ones(shuffled.size), 4)
    shuffled_sigma, _ = branching_parameter(
        surrogate["n_ancestors"].to_numpy(), surrogate["n_descendants"].to_numpy(), 60
    )
    _, _, shuffled_vuong = fit_mle(surrogate["size"].to_numpy(), 1, 60)

    check("shuffle collapses sigma", shuffled_sigma < 0.2, f"{shuffled_sigma:.4f}")
    check("shuffle destroys the power law", shuffled_vuong < -2, f"z={shuffled_vuong:.1f}")
    check("shuffle shrinks the largest avalanche",
          surrogate["size"].max() < found["size"].max() / 10,
          f"{surrogate['size'].max()} vs {found['size'].max()}")


def main():
    check_samplers()
    check_refractory()
    check_detection()
    check_branching()
    check_fitting()
    check_homogeneous_prediction()
    check_critical_control()

    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {', '.join(FAILURES)}")
        raise SystemExit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
