"""Generate an IID n x n random field over X iterations and write it to CSV.

The generator deliberately contains no spatial or temporal dependence: every
value is drawn independently. Any structure found later by analyze.py is
therefore feigned by construction.
"""

import argparse
import os

import numpy as np
import pandas as pd


def sample_powerlaw(rng, size, a, xmin, xmax):
    """Inverse-CDF sampling from the truncated PDF p(x) proportional to x**a."""
    u = rng.random(size)
    if np.isclose(a, -1.0):
        # CDF of 1/x integrates to a log, so the inverse is exponential.
        return xmin * (xmax / xmin) ** u
    if np.isclose(a, -2.0):
        return 1.0 / (1.0 / xmin - u * (1.0 / xmin - 1.0 / xmax))
    exponent = a + 1.0
    return ((xmax**exponent - xmin**exponent) * u + xmin**exponent) ** (1.0 / exponent)


def sample_gaussian(rng, size, mu, sigma):
    return np.clip(rng.normal(mu, sigma, size), 0.0, 100.0)


def simulate(args):
    rng = np.random.default_rng(args.seed)
    shape = (args.x, args.n, args.n)

    if args.mode == "powerlaw":
        if args.xmin <= 0:
            raise ValueError("--xmin must be positive for a power-law distribution")
        if args.xmax <= args.xmin:
            raise ValueError("--xmax must be greater than --xmin")
        values = sample_powerlaw(rng, shape, args.a, args.xmin, args.xmax)
    else:
        values = sample_gaussian(rng, shape, args.mu, args.sigma)

    iteration, row, col = np.meshgrid(
        np.arange(args.x), np.arange(args.n), np.arange(args.n), indexing="ij"
    )
    return pd.DataFrame(
        {
            "iteration": iteration.ravel(),
            "row": row.ravel(),
            "col": col.ravel(),
            "value": values.ravel(),
        }
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["powerlaw", "gaussian"], default="powerlaw")
    parser.add_argument("--n", type=int, default=16, help="grid side length")
    parser.add_argument("--x", type=int, default=200, help="number of iterations")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="data/timeseries.csv")
    parser.add_argument("--a", type=float, default=-2.0)
    parser.add_argument("--xmin", type=float, default=1.0)
    parser.add_argument("--xmax", type=float, default=100.0)
    parser.add_argument("--mu", type=float, default=50.0)
    parser.add_argument("--sigma", type=float, default=15.0)
    return parser.parse_args()


def main():
    args = parse_args()
    frame = simulate(args)

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    frame.to_csv(args.out, index=False, float_format="%.6f")

    print(f"mode={args.mode} n={args.n} iterations={args.x} seed={args.seed}")
    print(f"wrote {len(frame)} rows to {args.out}")
    print(
        f"value stats: min={frame['value'].min():.4f} "
        f"median={frame['value'].median():.4f} max={frame['value'].max():.4f}"
    )


if __name__ == "__main__":
    main()
