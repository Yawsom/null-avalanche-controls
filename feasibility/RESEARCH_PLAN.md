# Predicting avalanche statistics across bin widths

Prospective study design, 5 September 2026. This follows the completed feasibility
audit; the experiments below have not yet been run. It refines the research
question rather than replacing or retrospectively changing the original protocol.

## Question and potential contribution

How much of the bin-width dependence of measured avalanche statistics can be
predicted from event rate, shared fluctuations, refractoriness and observation
alone, and when does adding propagation improve predictions on held-out data?

The relationship between binning and apparent criticality is established.
[Priesemann & Shriki (2018)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006081)
derive avalanche distributions for homogeneous and slowly varying Poisson drive.
[Pausch (2022)](https://doi.org/10.1103/PhysRevResearch.4.023212) treats continuous-time
branching with immigration, binning effects and inference without time bins.
Our candidate contribution is a calibrated, finite-recording comparison of what
whole bin-width curves can identify under matched observation conditions, with
held-out empirical validation. Novelty and successful discrimination remain open.

An exact crossing of sigma=1 and tau=1.5 is not a success criterion. Failure to
fit an exact power law does not erase a reproducible bin-width relationship.

## 1. Derive and verify the baseline

Start with independent homogeneous Poisson population counts. For rate r and
width Delta, let mu=r*Delta, p=exp(-mu) and q=1-p. An avalanche is a maximal run
of occupied bins, and size counts events, not distinct electrodes.

- The expected repository single-ancestor estimate is mu, conditional on having
  eligible starts. Its selection is a blank bin followed by exactly one event;
  independence leaves the next-bin expectation unchanged.
- Duration L in bins satisfies P(L=l)=p*q^(l-1), l>=1, and E[L]=1/p.
- The probability-generating function of an occupied bin's event count is
  g(z)=(exp(mu*z)-1)/(exp(mu)-1). Consequently avalanche size has generating
  function F(z)=p*g(z)/(1-q*g(z)). This follows by summing over the geometric
  number of independent occupied bins. E[S]=mu/(p*q).

These equations give an exact null distribution and a check on implementation.
They assume no dead time, stationarity, independent increments and negligible
record-boundary censoring. Distinct-electrode sizes require a different count
law; refractory or shared-drive models require additional calculations.

Compute the expected fitted finite-support exponent from that distribution,
using exactly the estimator and support applied to simulations. For an interior
MLE, its population limit matches the conditional mean of log(size) to that of
the fitted power law. It is an estimator-dependent summary even when the true
distribution is not a power law. Validate analytic predictions against fresh
simulations and account for recordings with no eligible branching estimates.

Then extend to specified shared-drive models analytically where tractable and
with numerical predictions otherwise. Treat the existing results as exploratory
checks, not independent validation of choices they informed.

## 2. Compare mechanisms under matched conditions

Use three levels of model:

1. Independent activity, establishing rate and occupancy effects.
2. Shared time-varying drive with per-channel refractoriness, without propagation.
3. The same drive and observation model with explicit propagation added.

Use continuous-time events and distributed transmission delays to avoid building
a preferred 4 ms clock into the propagating model. Include smooth drive as well
as abrupt bursts, and vary drive strength and timescale independently of
propagation where feasible. Match observed rates, channel counts, recording
lengths and observation rules; report residual mismatches instead of silently
attributing their consequences to propagation.

Track causal ancestry and realized offspring in simulations. A nominal branching
parameter of one is insufficient to label a finite refractory model critical.
Subcritical propagation is an essential control. Apply the same subsampling and
event-count versus active-channel analyses to every family.

Use independent training and evaluation seeds. Choose nuisance ranges and sample
sizes using pilot precision and power estimates, then freeze the evaluation
protocol. Report detection and false-positive frequencies across those ranges,
including regimes where the models cannot be distinguished.

## 3. Predict the whole curve

Fit model parameters on training recordings using specified event-level
information: channel rates, refractory/ISI structure and population fluctuation
statistics at declared timescales. Fit no avalanche exponent or branching curve.
Specify which correlations enter fitting: predicting summaries of those same
correlations is not independent evidence of a new mechanism.

Simulate independent recordings from fitted models. Compare their predictions
with held-out curves of the single-ancestor branching estimate, size and duration
distributions, and finite-support MLE exponents. Flag poor power-law fit; retain
the underlying distributions and do not present such exponents as scaling laws.
Use alternative branching estimators as declared sensitivity analyses, keeping
their definitions separate.

Show widths in milliseconds and normalized by rate (r*Delta), with global and
conditional IEI markers. Test how much normalization explains and what residual
dependence remains on drive timescale and propagation. Rate normalization is not
expected to remove effects of changing refractory or drive timescale ratios.

Freeze widths, supports, bin-origin sensitivity and a joint discrepancy before
evaluation. Calibrate that discrepancy and simultaneous predictive bands using
whole simulated recordings, including parameter uncertainty. Different widths
from one recording are dependent observations, not independent replications.
Where likelihoods are tractable, add held-out event likelihood comparisons;
otherwise explicitly call the assessment a predictive summary comparison.

## 4. Validate on neural recordings

Select an accessible dataset with event timestamps, channel identities, documented
preprocessing and multiple independent preparations or sessions. No dataset has
yet been selected or analyzed for this phase. Begin with one modality; spike
results do not directly reanalyze LFP avalanches without an observation model.

Split by preparation/session where possible. For within-session prediction, use
separated time blocks and distinguish it from generalization across preparations.
Estimate drive dynamics and all tuning choices from training data. Reusing the
held-out population rate envelope is a conditional reconstruction and must be
reported separately from an independent held-out prediction.

Assess curve coverage and the incremental predictive value of propagation across
independent recordings, with uncertainty. Neural correlations make a structured
null necessary; their presence alone does not specify which mechanism generated
them. A single flexible null or a single poorly fitting alternative cannot settle
that question.

## Interpretation and milestones

| Result | Supported interpretation |
| --- | --- |
| Shared-drive model predicts held-out curves | These statistics are compatible with that nonpropagating explanation; they do not by themselves identify propagation. |
| Adding propagation reliably improves predictions | The curves contain information beyond the tested shared-drive family, under the stated modeling assumptions. This is not universal causal identification or proof of criticality. |
| Both families predict similarly | Quantify the parameter and recording-length regimes in which bin-width curves cannot distinguish them. |
| Both fail | Improve the drive or observation models; model failure is not evidence specifically for criticality. |

First milestone: verified analytic baseline and a matched synthetic comparison
showing either useful discrimination or a quantified region of indistinguishability.
Only then scale empirical validation. The completed feasibility code supplies
detectors, fitting checks and exploratory sweeps; it does not yet supply matched
model fitting or held-out empirical predictions.

The intended paper figures are: analytic versus simulated baseline curves;
matched model trajectories and discrimination across nuisance parameters; and
held-out neural trajectories with predictive bands. A credible manuscript needs
a useful result beyond the existing literature, not merely more parameter sweeps.
