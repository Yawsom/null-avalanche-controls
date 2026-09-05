# Literature and replication audit

Primary sources checked on 5 September 2026. This is a focused novelty screen,
not a systematic review. The interpretations in the last column are our
assessment, not the cited authors' claims about this repository.

| Source | Relevant result | Consequence for this project |
| --- | --- | --- |
| [Beggs & Plenz, 2003, Neuronal Avalanches in Neocortical Circuits](https://pmc.ncbi.nlm.nih.gov/articles/PMC6741045/) | Experimental avalanche statistics, branching estimates, spatial rescaling and pharmacological manipulation. | A simulation of selected temporal statistics cannot adjudicate the whole empirical argument. |
| [Touboul & Destexhe, 2010, Can Power-Law Scaling and Neuronal Avalanches Arise from Stochastic Dynamics?](https://pmc.ncbi.nlm.nih.gov/articles/PMC2820096/) | Stochastic surrogates can appear power-law-like under graphical analysis; stronger statistical tests distinguish some cases. | The regression-versus-statistical-fit caution is established. |
| [Klaus, Yu & Plenz, 2011, Statistical Analyses Support Power Law Distributions Found in Neuronal Avalanches](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019779) | Reanalyzes experimental data, including the 2003 culture data, using finite-size scaling, likelihood fits, and alternative distributions. | A critique based only on the lack of modern model comparison in 2003 misses a direct subsequent response. We have not reanalyzed those data. |
| [Touboul & Destexhe, 2017, Power-law statistics and universal scaling in the absence of criticality](https://arxiv.org/abs/1503.08033) | Noncritical systems and stochastic surrogates can exhibit scaling signatures. | Nonpropagating or noncritical counterexamples are not a new broad premise. |
| [Priesemann & Shriki, 2018, Can a time varying external drive give rise to apparent criticality in neural systems?](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006081) | Analytically treats homogeneous and time-varying Poisson activity. Shared-drive activity yields approximate avalanche power laws; shuffling removes them. Discusses additional statistics and limits of observational inference. | This is the closest overlap: mechanism, homogeneous control, shuffled comparison, and the interpretation problem are already covered. A different burst intensity distribution alone is a weak novelty claim. |
| [Pausch, 2022, From neuronal spikes to avalanches: Effects and circumvention of time binning](https://doi.org/10.1103/PhysRevResearch.4.023212) | Develops a continuous-time branching-with-immigration model, derives spike statistics, studies binning and subsampling, and compares with experimental data. | Binning effects and model-based inference without bins are already explicit research targets. Our prospective contribution must establish additional utility under realistic nuisance and observation conditions. |
| [Morrell, Nemenman & Sederberg, 2024, Neural criticality from effective latent variables](https://elifesciences.org/articles/89337) | Noninteracting units driven by latent variables exhibit size/duration scaling and exponent relationships across parameter ranges. | A contemporary benchmark must address more than sigma and a size exponent. |
| [Amarasingham et al., 2012, Conditional modeling and the jitter method of spike resampling](https://pubmed.ncbi.nlm.nih.gov/22031767/) | Interval jitter conditions on coarse counts and tests fine timing assumptions. Its response depends on the timescale of shared rate fluctuations. | Interval jitter is an existing method, and rejection is not evidence specifically for direct propagation. Our burst-boundary failures are compatible with this established limitation. |
| [Platkiewicz, Stark & Amarasingham, 2017, Spike-Centered Jitter Can Mistake Temporal Structure](https://pmc.ncbi.nlm.nih.gov/articles/PMC5955204/) | Gives false-positive examples for spike-centered jitter and distinguishes it from interval jitter. | Replacing one timing surrogate with another cannot be marketed as a new general correction without calibration. |
| [Clauset, Shalizi & Newman, Power-law distributions in empirical data](https://arxiv.org/abs/0706.1062) | Separates estimation, absolute fit testing and comparison to alternatives. | Our finite-support bootstrap is a conditional adaptation of established methodology, not a new statistical test. |

## Paper-to-code audit

The following source summary is intentionally limited to the main discrepancies.
In the [2003 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6741045/), events come
from filtered LFP maxima, with a 20 ms exclusion window; the repo directly
generates events and amplitudes. The paper derives the IEI cutoff from measured
cross-correlations (150–200 ms in cultures), whereas this repo fixes it at
200 ms. Both use consecutive occupied bins and an initial-bin branching
estimate. The paper's jitter experiment perturbs event timing; the repo also
permutes amplitudes. Its 70 hours combine seven cultures, whereas a single
70-hour simulation is one realization. The paper includes spatial rescaling and
pharmacological evidence that this repo does not reproduce. These differences
justify calling this a partial methodological reimplementation, not a faithful
replication or a refutation of the experimental paper.

## Additional repository-specific observations

These are code/math findings, not claims attributed to the papers above:

- `bp_simulate.generate_critical` bars reuse for the whole cascade. This is an
  additional finite-resource assumption, not a consequence of a 20 ms window.
- A finite electrode array alone does not cap total event count at 60 if sites
  can reactivate. The legacy model imposes that cap by construction.
- The label "homogeneous Bernoulli" is imprecise: generation uses Poisson counts
  with uniform integer times followed by dead-time thinning.
- Shared bursts can overlap. They are additive drive episodes, not a strict
  alternating on/off telegraph process. Uniform drive inside each episode plus
  abrupt edges is an important assumption for the jitter tests.
- For independent Poisson bin counts of mean lambda, occupancy is
  q=1-exp(-lambda). Conditional on an occupied start, a run's duration is
  geometric, P(L=l)=(1-q)q^(l-1), with mean exp(lambda). There is no singularity
  at lambda=1. Refractoriness and rate changes invalidate the simple IID premise.
- For IID bin counts, the single-ancestor estimator has expectation E[N_next],
  because N_next is independent of the previous blank bin and single-ancestor
  selection. This identity does not extend to arbitrary tissue dynamics or
  force unity with a conditional/truncated IEI definition.
- The amplitude normalization uses the smallest observed avalanche instead of
  the stored event threshold. The new experiments avoid amplitude inference.
- Count-based analysis can count repeated activations within wide bins, while
  the frame definition refers to active electrodes. Interpret bins beyond the
  refractory window as additional implementation sensitivity.
- Empty input and zero bin width are not safely handled by the original
  detector. The feasibility wrapper validates inputs; the original scripts are
  retained unchanged for a traceable baseline.

## Interpretation of novelty

The strongest potential contribution would be a validated benchmark of concrete
analysis decisions at finite recording lengths, tied to current workflows and
with demonstrable incremental utility. Neither re-showing shared-drive scaling
nor documenting that full shuffling destroys it establishes that contribution.
No claim of exhaustive absence of novelty is made.
