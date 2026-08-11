This repository holds two separate experiments that share a theme. Both ask
whether an analysis method reports structure in data that has none.

1. **[Power-Law Time Series Avalanche Experiment](#power-law-time-series-avalanche-experiment)**
   (`simulate.py`, `analyze.py`) — a similarity criterion applied to IID data.
2. **[Beggs-Plenz Replication on Null Data](#beggs-plenz-replication-on-null-data)**
   (`bp_simulate.py`, `bp_analyze.py`, `bp_validate.py`, `bp_compare.py`) — the
   neuronal avalanche method of Beggs & Plenz 2003 applied to non-propagating data.

---

# Power-Law Time Series Avalanche Experiment

Does a completely independent random process still look correlated when you
apply a neuroscience-style similarity criterion across neighbouring cells and
consecutive time bins?

The generator here has **no spatial or temporal dependence at all**. Every
value in the `n x n` grid, at every iteration, is drawn independently. So any
"propagation" the analysis reports is false by construction. What the
experiment measures is how strongly the choice of marginal distribution
(truncated power law vs Gaussian) inflates the rate of these feigned
correlations, and whether the false links assemble into avalanche-like events.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Simulate

```bash
python simulate.py --mode powerlaw --n 16 --x 200 --seed 42 --out data/timeseries_powerlaw.csv
python simulate.py --mode gaussian --n 16 --x 200 --seed 42 --out data/timeseries_gaussian.csv
```

Power-law mode samples the truncated PDF `p(x) ∝ x^a` on `[xmin, xmax]`
(defaults `a = -2`, `xmin = 1`, `xmax = 100`) by inverse-CDF transform of a
uniform draw. Gaussian mode samples `N(mu, sigma)` (defaults `50`, `15`) and
clips to `[0, 100]`.

Output is long form, with every frame retained:

```text
iteration,row,col,value
0,0,0,12.400000
```

### Simulation flags

| Flag | Default | Meaning |
| --- | --- | --- |
| `--mode` | `powerlaw` | `powerlaw` or `gaussian` |
| `--n` | `16` | grid side length |
| `--x` | `200` | number of iterations |
| `--seed` | `42` | RNG seed |
| `--out` | `data/timeseries.csv` | output CSV path |
| `--a`, `--xmin`, `--xmax` | `-2`, `1`, `100` | power-law parameters |
| `--mu`, `--sigma` | `50`, `15` | Gaussian parameters |

## Analyse

```bash
python analyze.py --csv data/timeseries_powerlaw.csv --z 1 --plot-dir plots/powerlaw/
python analyze.py --csv data/timeseries_gaussian.csv --z 1 --plot-dir plots/gaussian/
```

`z = 1` is illustrative only. The whole point is to repeat the analysis across
several tolerances, since `P(|X - Y| <= z)` behaves very differently for the
two marginals.

### What counts as a correlation

A **feigned correlation** links a source site `(r, c, t)` to a descendant
`(r', c', t+1)` when both hold:

1. The descendant is in the Moore neighbourhood of the source,
   `max(|r' - r|, |c' - c|) <= 1`. That includes staying in the same cell, so
   a signal may hold still, move orthogonally, or move diagonally.
2. `|v(r,c,t) - v(r',c',t+1)| <= z`.

**There are no within-frame links.** Two cells observed in the same iteration
are never connected, no matter how similar their values or how close they sit.
If every cell at iteration `t` held an identical value, that alone would
produce zero edges. The reasoning is that a genuinely propagating signal needs
one time bin to move between cells.

### Avalanches

Each participating space-time site `(r, c, t)` is a graph node and each
qualifying link is an edge. Edges point forward in time, but membership is
resolved by treating them as undirected and taking connected components
(union-find). An avalanche is one such component containing at least one edge.

A front may drift, branch, merge, stay put, or return to a cell it already
visited. Revisiting does not split the event: `A_t -> B_(t+1) -> A_(t+2)` is a
single avalanche.

Per-avalanche metrics written to `<plot-dir>/avalanches.csv`:

| Metric | Definition |
| --- | --- |
| `size` | number of participating space-time sites; the same cell at two different iterations counts twice |
| `duration` | `t_max - t_min + 1` |
| `spatial_extent` | number of unique `(r, c)` cells recruited |
| `start_iteration`, `end_iteration` | `t_min`, `t_max` |
| `edge_count` | number of feigned-correlation edges in the component |

For `A_t -> B_(t+1) -> A_(t+2)`: `size = 3`, `duration = 3`,
`spatial_extent = 2`.

### Edge-rate diagnostic

Before any avalanche grouping, `analyze.py` reports

```text
R(z) = (candidate temporal pairs with |Δv| <= z) / (all candidate temporal pairs)
```

where candidates are every in-bounds Moore pair between consecutive frames.
This is the mechanism knob: it exposes how much of the apparent structure comes
purely from the marginal distribution. A sweep over `z` is written to
`<plot-dir>/edge_rate_vs_z.csv` and plotted alongside, controlled by
`--z-sweep` (comma-separated, default `0.01,0.03,0.1,0.3,1,3,10,30`).

## Outputs

```text
plots/<mode>/
├── avalanches.csv
├── avalanches_vs_size.png
├── avalanches_vs_duration.png
├── avalanches_vs_spatial_extent.png
├── edge_rate_vs_z.csv
└── edge_rate_vs_z.png
```

Distribution plots use log-log axes over raw counts. A straight-ish line on
log-log axes is **not** evidence of a power law; fitting and comparing against
alternatives would be a separate exercise.

## Interpreting the comparison

Both modes go through the identical analysis pipeline, so the only manipulated
variable is the IID marginal. The truncated power law with `a = -2` piles most
of its probability near `xmin`, which makes two independent draws land within
`z` of each other far more often than under a Gaussian spread across the same
range. The chain being tested:

```text
IID marginal -> P(|X - Y| <= z) -> feigned links -> connected space-time
events -> observed avalanche size and duration distributions
```

---

# Beggs-Plenz Replication on Null Data

Beggs & Plenz (2003) report that cortical activity organises into "neuronal
avalanches": a branching parameter `sigma` of 1.04 and a size distribution
following a power law with exponent `tau = 3/2`, together taken as evidence
that cortex operates at criticality.

Their avalanche definition has one free parameter, the bin width `dt`, which
they set to `IEI_avg`, the mean interval between events anywhere on the array.
This pipeline reimplements their method faithfully and runs it on data built to
contain **no propagation whatsoever**, sweeping `dt` to separate results that
reflect the data from results that reflect the binning choice.

## The headline result

A bursty null with zero electrode-to-electrode propagation reproduces the
published signature exactly:

| Condition | `dt` | `sigma` | `tau` | Beats an exponential? |
| --- | --- | --- | --- | --- |
| Beggs & Plenz 2003 | `IEI_avg` ≈ 4 ms | 1.04 ± 0.19 | 1.5 | not tested |
| **Bursty null, mixed intensity** | **10 ms** | **0.986** | **1.526** | **yes, z = +34.5** |
| Critical branching (positive control) | 4 ms | 1.001 | 1.478 | yes, z = +152.7 |
| Homogeneous null | any | reaches 1 at 64 ms | never 1.5 together | **no, at any `dt`** |
| Bursty null, uniform intensity | any | reaches 1 at 5 ms | never 1.5 together | **no, at any `dt`** |

In that bursty null every electrode is conditionally independent given a shared
global drive. Nothing propagates. The `(sigma, tau) = (1, 3/2)` signature still
appears, and the size distribution genuinely beats an exponential.

![phase comparison](plots/bp_comparison/phase_comparison.png)

Position in this plane is not the evidence: every trajectory sweeps past the
critical point at some bin width. The filled markers are what matters, showing
where the size distribution actually outperforms an exponential.

## What this does and does not say about the paper

Three findings cut in the paper's favour:

- **A homogeneous Poisson null never produces a power law.** Its Vuong statistic
  is negative at every bin width tested, from 1 to 128 ms. Independence alone is
  not enough to manufacture `tau = 3/2`.
- **Bursting alone is not enough either.** A telegraph drive with uniform burst
  intensity, matched to their `IEI_avg` of about 4.2 ms, still fails at every bin
  width. Reproducing the exponent needed a *broad spread of burst intensities*,
  not merely synchrony.
- **The paper's own bin rule does not select the artifact.** At `dt = IEI_avg`
  the mixed-intensity null gives `sigma = 0.605` and `tau = 2.02`, nowhere near
  the signature. Reaching it required `dt = 10 ms`, roughly 2.2x `IEI_avg`. Their
  prescription is not free to be tuned, and on this dataset it would not have
  produced the false positive.

Two findings cut against it:

- **`sigma` is close to worthless on its own.** For independent data the
  branching parameter is exactly the mean number of events per bin, so
  `sigma = rate x dt`. Binning at the mean inter-event interval forces
  `sigma = 1` by arithmetic. `bp_validate.py` confirms this to within 1 percent
  at four bin widths and shows independent data hitting `sigma = 0.9956` at
  `dt = 62 ms` while remaining decisively non-power-law. Their measured 1.04 is
  what the binning convention guarantees, whatever the tissue is doing.
- **Their jitter control cannot detect a false power law.** Jittering the
  positive control by 4 ms drops `sigma` from 1.00 to 0.61, yet the size
  distribution still beats an exponential at z = +196; at 80 ms jitter it is
  still z = +314. The paper applies jitter only to `sigma`, and our results show
  that is the only quantity it tests. A full time shuffle does destroy the
  signature (z = +153 to -180), so it is the far stronger control.

The net position: the branching parameter is an artifact of the binning
convention, but `tau = 3/2` is a real constraint that most nulls fail. It is
not, however, unique to propagating systems, and the paper does not report a
model comparison that would separate a power law from the alternatives.

## Two methodological cautions

**A straight line on log-log axes is not a power law.** Exponentially
distributed samples fed to log-log regression return "exponents" of 1.75 and
3.58 with R² above 0.81, while the Vuong test correctly rejects them. The MLE is
no safer used alone: on homogeneous null data it returns `tau = 1.398`, which
reads like 3/2, while the model comparison rejects the power law at z = -405.
Only the comparison is informative.

**Fits must stop short of the array size.** A bounded array puts a spike of
probability at exactly its own limit, since every cascade that would have grown
larger is truncated to land there. In the positive control, size 60 carries 1541
counts against about 10 for sizes 53-59. Fitting across that spike returns
`tau = 1.41`; fitting below it returns 1.478 with R² = 0.9998. Both are reported,
as `tau_regression` and `tau_regression_to_cutoff`.

## Generate

```bash
python bp_simulate.py --mode homogeneous --hours 70 --out data/bp_homogeneous.npz
python bp_simulate.py --mode bursty      --hours 70 --out data/bp_bursty.npz
python bp_simulate.py --mode bursty --burst-intensity-sigma 1.0 --hours 70 \
    --out data/bp_bursty_het.npz
python bp_simulate.py --mode critical    --hours 70 --out data/bp_critical.npz
```

Each run produces about 4 million events over 70 simulated hours, written as a
compressed `.npz` of roughly 37 MB holding `electrode`, `time_ms`, `amplitude`
and a JSON metadata blob.

### Generator modes

| Mode | Propagation | Description |
| --- | --- | --- |
| `homogeneous` | none | every electrode fires independently at a constant rate |
| `bursty` | none | a shared telegraph drive gates conditionally independent electrodes |
| `critical` | yes | Galton-Watson cascades; the positive control |

`--burst-intensity-sigma` is the knob that matters for the bursty null. A run
inside a burst of occupancy `lambda` is geometric with scale `1/(1 - lambda)`
below one and effectively unbounded above one, so spreading `lambda` across that
boundary mixes wildly different scales. The spread is lognormal, deliberately
chosen as broad but not itself a power law, so nothing was injected by hand.
`--burst-duration-sigma` spreads burst lengths instead and barely matters, since
an avalanche ends at the first empty bin regardless of how long the burst runs.

### Parameters matched to the paper

| Quantity | Value |
| --- | --- |
| Electrodes | 60 (8x8 minus corners) |
| Time resolution | 1 ms |
| Refractory period | 20 ms per electrode |
| Array-wide rate | 58,000 events/hour |
| Burst duty cycle | 6.8 percent, ~150 ms bursts, giving `IEI_avg` ≈ 4.5 ms against their 4.2 ms |

Amplitudes are the supra-threshold values, drawn exactly from the tail of the
chosen marginal. For a power law the tail is another power law with the same
exponent, so amplitudes span 3,710 to 999,890, matching the several orders of
magnitude the paper describes. A Gaussian marginal at the same event rate yields
a threshold of 3.46 SD and amplitudes spanning only 3.46 to 6.09, which is why
`--amp-dist` changes the amplitude-weighted result so sharply.

## Analyse

```bash
python bp_analyze.py --events data/bp_bursty_het.npz --plot-dir plots/bp_bursty_het/
python bp_analyze.py --events data/bp_critical.npz --native-bin-ms 4 \
    --plot-dir plots/bp_critical/
python bp_analyze.py --events data/bp_critical.npz --control shuffle \
    --plot-dir plots/bp_critical_shuffle/

python bp_compare.py     # overlays every sweep onto one comparison
python bp_validate.py    # 50 correctness checks; run before trusting anything
```

### What counts as an avalanche

Their definition is purely temporal, and spatial adjacency plays no part; the
paper itself reports a contiguity index of only 39 percent, meaning activity
usually skips its nearest neighbours. Events are binned at width `dt`, and an
avalanche is a maximal run of consecutively occupied bins bracketed by empty
ones. A single occupied bin between two empty ones is a valid avalanche.

| Metric | Definition |
| --- | --- |
| `size` | total electrode activations in the run, counting repeats, which is why sizes can exceed 60 |
| `duration` | number of bins spanned |
| `amplitude_size` | summed amplitude, their second size definition |
| `sigma_single` | mean descendants per ancestor across the first two bins, restricted to single-ancestor avalanches |
| `sigma_multi` | ancestor-weighted variant carrying their availability correction `(nmax - 1) / (nmax - na)` |

### Analysis flags

| Flag | Default | Meaning |
| --- | --- | --- |
| `--events` | required | `.npz` from `bp_simulate.py` |
| `--plot-dir` | required | output directory |
| `--bins-ms` | `1,...,128` | bin widths to sweep |
| `--tmax-ms` | `200` | truncation for the conditional `IEI_avg` |
| `--native-bin-ms` | `0` | summarise at this width instead of the `IEI_avg` one |
| `--fit-max` | `0` | largest size in the fit; 0 uses half the electrode count |
| `--control` | `none` | `none`, `jitter` or `shuffle` |
| `--jitter-ms` | `4` | jitter magnitude |

`IEI_avg` is recomputed per dataset as the conditional mean under `--tmax-ms`,
replicating their procedure rather than hardcoding 4.2 ms. Controls re-enter the
identical binning path and have `IEI_avg` recomputed on the surrogate.

## Outputs

```text
plots/bp_<condition>/
├── bin_sweep.csv                     # one row per dt: sigma, tau, Vuong, counts
├── bin_sweep.png
├── avalanches_at_iei.csv
├── size_distribution.png
├── duration_distribution.png
├── amplitude_size_distribution.png
└── phase_trajectory.png

plots/bp_comparison/
├── phase_comparison.png              # all conditions, Figure 7D style
├── vuong_vs_binwidth.png
└── summary.csv
```

## A caveat on the generators

The bursty modes lose events to the refractory period, since bursts concentrate
activity: the uniform-intensity null realises 53,804 events/hour against the
58,000 target, and the mixed-intensity null 49,549. Rates are reported on every
run. The deficit is not corrected because `IEI_avg` is recomputed per dataset, so
each condition is analysed at its own appropriate bin width, and the spread is
far smaller than the 10,000 to 240,000 per hour range across their own cultures.

The positive control bars an electrode from being recruited twice within one
cascade. This follows from the 20 ms refractory period, as avalanches are much
shorter than that, and it supplies the cutoff at the array size that the paper
measures. It matters for a second reason: at `sigma` exactly 1 an unrestricted
critical process wanders unboundedly before dying, and on a 60-electrode array
it settles into a self-sustaining state producing single cascades of over
170,000 events. Exhausting the array is what actually terminates it.
