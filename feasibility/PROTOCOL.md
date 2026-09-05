# Publication feasibility protocol — 5 September 2026

This is a bounded computational feasibility study, not a preregistered study or
a completed paper. It was written before running the multi-seed experiments,
after the exploratory seed-42 results documented in the conversation. Those
results informed the choices below. Original experiment scripts remain intact.

## Questions and decision

1. Does the reported mixed-drive null reproduce the joint signature under the
   repository's IEI rule, rather than only after bin selection?
2. Does the result survive consistent exponent estimation, absolute fit checks,
   changes in fitted support, and independent simulation seeds?
3. Does loss of structure after shuffling distinguish shared drive from
   propagation? Can an existing local timing randomization improve separation?
4. Does this add enough to the closest primary literature to justify a paper?

Proceed to a manuscript only if the study identifies a robust, distinct failure
mode with a useful consequence for practice. Otherwise recommend a reproducible
technical note or stop. A null result or a failed proposed diagnostic is a valid
outcome. Do not tune parameters to manufacture a positive feasibility decision.

## Experiment grid

- Main replications: seeds 100–119, 0.25 and 3 simulated hours, 60 electrodes.
  Original homogeneous, uniform-burst, mixed-burst (log-intensity SD 1), and
  cascade-capped branching generators.
- Positive-control sensitivity: a finite array branching model with an actual
  20 ms refractory window, allowing reuse after recovery; attempted offspring
  means 0.8, 1.0, 1.2. These are nominal parameters, not verified critical states.
  Apply the same durations and seeds as above. Record effective reproduction,
  event rate, cascade duration, and censoring/termination counts.
- Long-run confirmation: original mixed-burst model, seeds 42 and 100–104,
  70 simulated hours. This does not represent seven biological preparations.
- Main bin sweep: 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64, 96, 128 ms,
  plus the rounded conditional IEI. The IEI includes simultaneous events and
  conditions on gaps <=200 ms, exactly as in the original repository. Also
  inspect unrounded IEI and conditioning cutoffs 100 and 400 ms.
- Primary fitting support: integer sizes 1–30. Sensitivity: 1–15, 1–59, 3–30.
  These are conditional finite-support fits, not evidence for an asymptotic tail.
- A descriptive signature uses |sigma−1|<0.2 and |tau−1.5|<0.2, matching the
  original comparison tolerance. Report regression and MLE separately. Retain
  the original power-law-versus-exponential z>2 criterion only as a descriptive
  legacy screen. It is not a calibrated causal or criticality test.
- Report prescribed IEI, fixed 10 ms, and best bin from the stated sweep
  separately. Selecting a bin and testing on the same recording is exploratory.

## Statistical checks

- Fit discrete finite-support power laws and exponentials by likelihood. Add a
  lognormal-shaped discrete family on exactly the same support, using
  p(s) proportional to exp(−a log(s)−b log(s)^2), b>=0. It includes the power law
  as b=0. Report AIC differences; do not apply a nonnested Vuong test to this
  boundary/nested comparison. AIC is descriptive when samples are dependent.
- Absolute fit: discrete KS statistic with 199 parametric multinomial bootstrap
  replicates, re-estimating the exponent in each replicate. Use (1+exceedances)/200.
  This is an IID-model compatibility diagnostic; avalanche dependence can make
  its nominal p-values anti-conservative. Report independent-run variation as
  the primary replication evidence and add 1-second block-bootstrap checks on
  selected long-run/null and positive-control records.
- Include a known exact finite-support power-law sampler and an exact critical
  Poisson Galton–Watson total-progeny (Borel) sampler to check calibration and the
  distinction between asymptotic scaling and an exact power law on sizes 1–30.
- Record the selected sample size and fraction of all avalanches fitted.

## Surrogate experiment

On the 0.25-hour main recordings, retain the observed IEI bin and fixed 10 ms for
paired full-shuffle and ±4/±20/±80 ms jitter checks. Do not permute amplitudes;
the feasibility study uses event counts only. These surrogate trains need not
preserve refractoriness and that limitation must be reported.

Test a candidate short-lag diagnostic on 4 ms population counts: sum of products
of adjacent-bin counts. Compare it with 99 interval-jitter surrogates preserving
each electrode's event count in fixed 20, 40, and 100 ms windows. Use a one-sided
Monte Carlo rank p-value <=0.05. These surrogates assume timing exchangeability
inside windows, which is not guaranteed for burst boundaries or refractory
trains. Measure their actual false-positive frequency on the simulated nulls;
do not claim they test causality. Report all window sizes, not only the best.

## Reproducibility and limits

Save the run configuration, source hashes, environment versions, per-run CSVs,
summary tables, figures, tests, and an evidence-backed written decision under
`feasibility/`. Do not retain large raw event arrays; deterministic seeds and
generator code reproduce them. This audit does not use experimental neural data,
reproduce LFP preprocessing, or establish fidelity of every step to the 2003
paper. It cannot settle whether cortex is critical.
