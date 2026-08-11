"""Overlay the bin-width sweeps from every condition onto one comparison.

Reads the bin_sweep.csv files written by bp_analyze.py and produces the figure
that actually answers the question: for each dataset, does the trajectory
through (sigma, tau) space pass through the critical point the paper reports,
and does the size distribution survive a comparison against an exponential
when it gets there.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PAPER_SIGMA = 1.0
PAPER_TAU = 1.5

CONDITIONS = [
    ("homogeneous", "plots/bp_homogeneous/bin_sweep.csv", "#718096", "Homogeneous (null)"),
    ("bursty", "plots/bp_bursty/bin_sweep.csv", "#2b6cb0", "Bursty, uniform intensity (null)"),
    ("bursty_het", "plots/bp_bursty_het/bin_sweep.csv", "#805ad5", "Bursty, mixed intensity (null)"),
    ("critical", "plots/bp_critical/bin_sweep.csv", "#c53030", "Critical branching (positive control)"),
]


def signature_distance(frame):
    return np.hypot(frame["sigma_single"] - PAPER_SIGMA, frame["tau_regression"] - PAPER_TAU)


def load():
    loaded = []
    for key, path, colour, label in CONDITIONS:
        if not os.path.exists(path):
            print(f"skipping {key}: {path} not found")
            continue
        frame = pd.read_csv(path)
        frame["condition"] = key
        loaded.append((key, frame, colour, label))
    return loaded


def plot_phase(loaded, path):
    figure, axes = plt.subplots(figsize=(7.5, 5.5))

    # Every trajectory sweeps past the critical point, so position alone proves
    # nothing. Filled markers are the bin widths where the size distribution
    # actually beats an exponential; hollow ones only look the part.
    for index, (_, frame, colour, label) in enumerate(loaded):
        axes.plot(frame["sigma_single"], frame["tau_regression"], "-",
                  color=colour, linewidth=1.2, label=label, alpha=0.85)

        genuine = frame["vuong_z"] > 2
        axes.scatter(frame.loc[genuine, "sigma_single"], frame.loc[genuine, "tau_regression"],
                     s=34, color=colour, zorder=4)
        axes.scatter(frame.loc[~genuine, "sigma_single"], frame.loc[~genuine, "tau_regression"],
                     s=30, facecolors="white", edgecolors=colour, linewidths=1.1, zorder=4)

        best = frame.loc[signature_distance(frame).idxmin()]
        axes.annotate(f"{best['dt_ms']:.0f} ms", (best["sigma_single"], best["tau_regression"]),
                      fontsize=8, color=colour, fontweight="bold",
                      xytext=(8, 8 - 16 * (index % 2)), textcoords="offset points")

    axes.scatter([], [], s=34, color="#4a5568", label="filled: beats an exponential")
    axes.scatter([], [], s=30, facecolors="white", edgecolors="#4a5568",
                 label="hollow: exponential fits better")

    axes.scatter([PAPER_SIGMA], [PAPER_TAU], marker="*", s=340, color="#1a202c", zorder=6,
                 label="Beggs & Plenz (1.0, 1.5)")
    axes.axvline(PAPER_SIGMA, linewidth=0.6, alpha=0.35)
    axes.axhline(PAPER_TAU, linewidth=0.6, alpha=0.35)

    axes.set_xlim(0, 3)
    axes.set_ylim(0, 4)
    axes.set_xlabel("branching parameter sigma")
    axes.set_ylabel("size exponent tau (log-log regression)")
    axes.set_title("Each curve is one dataset swept over bin width\n"
                   "labels mark the bin width closest to the reported critical point",
                   fontsize=10)
    axes.legend(frameon=False, fontsize=8, loc="upper right")
    axes.grid(True, linewidth=0.3, alpha=0.4)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def plot_vuong(loaded, path):
    """Whether the size distribution actually beats an exponential, versus dt."""
    figure, axes = plt.subplots(figsize=(7.5, 4.5))

    for _, frame, colour, label in loaded:
        axes.plot(frame["dt_ms"], frame["vuong_z"], "-o", color=colour,
                  markersize=3.5, linewidth=1.2, label=label, alpha=0.85)

        # At very wide bins nearly everything merges into a few huge avalanches
        # and the fitted exponent goes negative. Such a fit beats an exponential
        # only because the data is increasing, which is not a power law in any
        # sense the paper means, so it is flagged rather than counted.
        degenerate = frame["tau_regression"] <= 0
        axes.scatter(frame.loc[degenerate, "dt_ms"], frame.loc[degenerate, "vuong_z"],
                     marker="x", s=55, color="#1a202c", linewidths=1.4, zorder=5)

    axes.scatter([], [], marker="x", s=55, color="#1a202c", linewidths=1.4,
                 label="degenerate fit (tau <= 0)")
    axes.axhline(0, color="#1a202c", linewidth=0.8)
    axes.axhspan(-2, 2, color="#a0aec0", alpha=0.25)
    axes.text(1.1, 0, " inconclusive", fontsize=7, va="center", color="#4a5568")
    axes.set_xscale("log")
    axes.set_yscale("symlog", linthresh=10)
    axes.set_xlabel("bin width dt (ms)")
    axes.set_ylabel("Vuong z: power law vs exponential")
    axes.set_title("Positive means the power law fits better than an exponential", fontsize=10)
    axes.legend(frameon=False, fontsize=8)
    axes.grid(True, which="both", linewidth=0.3, alpha=0.4)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def summarise(loaded):
    rows = []
    for key, frame, _, label in loaded:
        native = frame.loc[(frame["dt_ms"] - frame["iei_avg_ms"]).abs().idxmin()]
        closest = frame.loc[signature_distance(frame).idxmin()]
        rows.append({
            "condition": key,
            "description": label,
            "iei_avg_ms": round(float(frame["iei_avg_ms"].iloc[0]), 2),
            "at_iei_dt_ms": int(native["dt_ms"]),
            "at_iei_sigma": round(float(native["sigma_single"]), 3),
            "at_iei_tau": round(float(native["tau_regression"]), 3),
            "at_iei_vuong_z": round(float(native["vuong_z"]), 1),
            "best_dt_ms": int(closest["dt_ms"]),
            "best_sigma": round(float(closest["sigma_single"]), 3),
            "best_tau": round(float(closest["tau_regression"]), 3),
            "best_vuong_z": round(float(closest["vuong_z"]), 1),
            "reproduces_signature": bool(
                abs(closest["sigma_single"] - PAPER_SIGMA) < 0.2
                and abs(closest["tau_regression"] - PAPER_TAU) < 0.2
                and closest["vuong_z"] > 2
            ),
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="plots/bp_comparison")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    loaded = load()
    if not loaded:
        raise SystemExit("no sweep files found; run bp_analyze.py first")

    plot_phase(loaded, os.path.join(args.out_dir, "phase_comparison.png"))
    plot_vuong(loaded, os.path.join(args.out_dir, "vuong_vs_binwidth.png"))

    summary = summarise(loaded)
    summary.to_csv(os.path.join(args.out_dir, "summary.csv"), index=False)

    pd.set_option("display.width", 200)
    print(summary.to_string(index=False))
    print(f"\nwrote {args.out_dir}/")


if __name__ == "__main__":
    main()
