"""Statistical and simulation helpers for the bounded feasibility audit.

Finite-support goodness of fit is an IID compatibility diagnostic. Neither
its p-value nor a surrogate rank test is a test of biological causation.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from bp_analyze import detect_avalanches, branching_parameter, fit_regression
from bp_validate import build_args
from bp_simulate import GENERATORS, apply_refractory


def scalar_family(counts, statistic, lower=-5.0, upper=12.0):
    """MLE for p(x) proportional to exp(-theta * statistic(x)); batched."""
    counts = np.atleast_2d(np.asarray(counts, float))
    statistic = np.asarray(statistic, float)
    totals = counts.sum(axis=1)
    if np.any(totals <= 0):
        raise ValueError('positive sample size required')
    target = (counts * statistic).sum(axis=1) / totals
    lo, hi = np.full(totals.size, lower), np.full(totals.size, upper)
    for _ in range(48):
        mid = (lo + hi) / 2
        logp = -mid[:, None] * statistic
        logp -= logsumexp(logp, axis=1)[:, None]
        expectation = (np.exp(logp) * statistic).sum(axis=1)
        lo = np.where(expectation > target, mid, lo)
        hi = np.where(expectation > target, hi, mid)
    theta = (lo + hi) / 2
    logp = -theta[:, None] * statistic
    logp -= logsumexp(logp, axis=1)[:, None]
    return theta, logp


def size_counts(sizes, low=1, high=30):
    sizes = np.asarray(sizes)
    if low < 1 or high <= low or np.any(sizes != np.floor(sizes)):
        raise ValueError('integer sizes and 1 <= low < high required')
    selected = sizes[(sizes >= low) & (sizes <= high)].astype(int)
    return np.bincount(selected - low, minlength=high-low+1), selected


def fit_counts(counts, low=1):
    support = np.arange(low, low + len(counts), dtype=float)
    theta, lp = scalar_family(counts, np.log(support))
    decay, le = scalar_family(counts, support, lower=1e-6, upper=5.)
    lp, le = lp[0], le[0]
    n = counts.sum()
    diff = lp - le
    llr = float(counts @ diff)
    variance = float(counts @ (diff - llr/n)**2 / n)
    return dict(tau_mle=float(theta[0]), exp_rate=float(decay[0]),
                llr=llr, vuong_z=llr / np.sqrt(n * variance) if variance else np.nan,
                ks=float(np.max(np.abs(np.cumsum(counts/n) - np.cumsum(np.exp(lp))))),
                fit_n=int(n)), lp


def absolute_fit(counts, low, rng, bootstrap=199):
    fit, lp = fit_counts(counts, low)
    n = int(counts.sum())
    draws = rng.multinomial(n, np.exp(lp), size=bootstrap)
    support = np.arange(low, low+len(counts))
    _, boot_lp = scalar_family(draws, np.log(support))
    distances = np.max(np.abs(np.cumsum(draws/n, axis=1) -
                              np.cumsum(np.exp(boot_lp), axis=1)), axis=1)
    fit['ks_iid_p'] = float((1 + np.count_nonzero(distances >= fit['ks'])) / (bootstrap+1))

    # Nested lognormal-shaped family. b=0 is exactly the power-law boundary.
    logs = np.log(support)
    features = np.stack((logs, logs**2), axis=1)
    empirical = counts @ features / n
    def objective(params):
        z = -features @ params
        norm = logsumexp(z)
        probabilities = np.exp(z-norm)
        return float(norm + empirical @ params), empirical - probabilities @ features
    result = minimize(objective, [fit['tau_mle'], 0.], jac=True,
                      method='L-BFGS-B', bounds=[(-100., 100.), (0., 100.)],
                      options={'ftol': 1e-12, 'gtol': 1e-8, 'maxiter': 1000})
    if not result.success:
        # L-BFGS can report a failed line search at a near-stationary boundary.
        # Use a different constrained optimizer, never silently accept failure.
        result = minimize(objective, result.x, jac=True, method='SLSQP',
                          bounds=[(-100.,100.),(0.,100.)],
                          options={'ftol':1e-12,'maxiter':1000})
    if not result.success:
        raise RuntimeError(f'lognormal optimization failed: {result.message}')
    fit['lognormal_a'], fit['lognormal_b'] = map(float, result.x)
    fit['lognormal_bound_hit'] = bool(abs(result.x[0]) > 99.9 or result.x[1] > 99.9)
    # Positive favors the lognormal-shaped family; two parameters versus one.
    fit['aic_powerlaw_minus_lognormal'] = float(2*(-n*result.fun - counts@lp)-2)
    return fit


def measure(times, dt):
    times = np.asarray(times)
    if not np.isfinite(dt) or dt <= 0 or times.ndim != 1 or times.size == 0:
        raise ValueError('nonempty one-dimensional times and positive finite bin required')
    if not np.all(np.isfinite(times)) or np.any(times < 0):
        raise ValueError('nonnegative finite times required')
    return detect_avalanches(times, np.ones(times.size), dt)


def metrics(found, low=1, high=30):
    sizes = found['size'].to_numpy()
    counts, _ = size_counts(sizes, low, high)
    sigma, _ = branching_parameter(found.n_ancestors.to_numpy(),
                                    found.n_descendants.to_numpy(), 60)
    tau, r2 = fit_regression(sizes, low, high)
    fit, _ = fit_counts(counts, low) if counts.sum() >= 50 else ({}, None)
    return dict(sigma_single=sigma, tau_regression=tau, regression_r2=r2,
                n_avalanches=len(found), fit_fraction=float(counts.sum()/len(found)),
                mean_size=float(sizes.mean()), max_size=int(sizes.max()), **fit)


def finite_refractory(rng, hours, nominal, electrodes=60, refractory_ms=20):
    """Finite-population branching, with reuse after recovery and collisions.

    Each active unit attempts Poisson(nominal) offspring onto uniformly chosen
    array units. Occupied/refractory targets and simultaneous collisions reduce
    realized reproduction. Nominal=1 is not asserted to be this model's critical
    point. Cascades are seeded after silence; no within-cascade immigration.
    """
    total_ms = int(hours * 3_600_000)
    target = int(hours * 58_000)
    cascades, produced, parents, children, attempted, terminated = [], 0, 0, 0, 0, 0
    durations = []
    while produced < target:
        active = np.array([rng.integers(electrodes)])
        last = np.full(electrodes, -refractory_ms)
        last[active] = 0
        ids, offsets = [active], [np.zeros(1, dtype=int)]
        generation = 0
        while active.size and generation < 10_000:
            wanted = int(rng.poisson(nominal * active.size))
            candidates = np.unique(rng.integers(0, electrodes, wanted))
            now = (generation+1)*4
            new = candidates[now-last[candidates] >= refractory_ms]
            parents += active.size
            attempted += wanted
            children += new.size
            generation += 1
            if new.size:
                last[new] = now
                ids.append(new)
                offsets.append(np.full(new.size, generation, dtype=int))
            active = new
        terminated += int(active.size > 0)
        offsets, ids = np.concatenate(offsets), np.concatenate(ids)
        cascades.append((offsets, ids))
        produced += ids.size
        durations.append(int(offsets.max())+1)
    durations = np.array(durations)
    spare = max(5*len(cascades), total_ms//4-int(durations.sum()))
    gaps = 5 + rng.poisson(max(0., spare/len(cascades)-5.), len(cascades))
    starts = np.cumsum(gaps + np.r_[0, durations[:-1]])
    phase = rng.integers(0, 4, len(cascades))
    times = np.concatenate([(start+off)*4+p for start,p,(off,_) in zip(starts,phase,cascades)])
    ids = np.concatenate([ids for _,ids in cascades]).astype(np.int16)
    inside = times < total_ms
    # Phase offsets can shorten the gap by up to 3 ms; minimum silence is 5 bins
    # plus the previous duration, so cross-cascade separation remains >=21 ms.
    return ids[inside], times[inside], dict(
        nominal_offspring=nominal, actual_offspring=children/parents,
        attempted_offspring=attempted/parents, generation_limit_hits=terminated,
        generated_cascades=len(cascades), max_generated_duration_ms=int(durations.max()*4),
        generated_fraction_over_20ms=float(np.mean(durations*4 > 20)),
        clipped_events=int((~inside).sum()))


def generate(condition, hours, seed):
    rng = np.random.default_rng(seed)
    if condition.startswith('refractory_'):
        ids, times, metadata = finite_refractory(rng, hours, float(condition.split('_')[1]))
    else:
        mode = {'homogeneous':'homogeneous', 'uniform':'bursty',
                'mixed':'bursty', 'capped':'critical'}[condition]
        args = build_args(mode=mode, burst_intensity_sigma=1. if condition=='mixed' else 0.)
        ids, times = GENERATORS[mode](rng,args,int(hours*3_600_000),58_000/3_600_000/60)
        order, keep = apply_refractory(ids,times,20)
        ids,times = ids[order][keep],times[order][keep]
        order = np.argsort(times, kind='stable')
        ids,times = ids[order],times[order]
        metadata = {'dropped_events':int((~keep).sum())}
    span = max(int(hours*3_600_000), int(times.max())+1)
    metadata.update(events=int(times.size), span_ms=span,
                    realized_rate_per_hour=times.size*3_600_000/span)
    return ids,times,metadata


def surrogate(times, rng, span, kind, width=4):
    if kind == 'shuffle':
        moved = rng.integers(0,span,times.size)
    elif kind == 'jitter':
        moved = np.clip(times+rng.integers(-width,width+1,times.size),0,span-1)
    elif kind == 'interval':
        left = times//width*width
        moved = left + rng.integers(0, np.minimum(width,span-left))
    else:
        raise ValueError(kind)
    # Each event keeps its electrode association before sorting. Population
    # statistic uses only times; the per-electrode/window counts are preserved.
    return np.sort(moved)


def lag_product(times, dt=4):
    bins, counts = np.unique(times//dt,return_counts=True)
    adjacent = np.diff(bins)==1
    return int(np.sum(counts[:-1][adjacent]*counts[1:][adjacent]))


def interval_test(times, rng, span, width, reps=99):
    observed = lag_product(times)
    scores = np.array([lag_product(surrogate(times,rng,span,'interval',width))
                       for _ in range(reps)])
    return dict(observed=observed, surrogate_mean=float(scores.mean()),
                surrogate_sd=float(scores.std(ddof=1)),
                excess=float(observed-scores.mean()),
                rank_p=float((1+np.count_nonzero(scores>=observed))/(reps+1)))


def block_fit(found, dt, span, rng, low=1, high=30, block_ms=1000, reps=199):
    """Ordinary time-block bootstrap with centered, refitted CDF residuals.

    Whole avalanches are assigned by their start time. Fixed disjoint blocks
    preserve local dependence approximately; long-range dependence and
    nonstationarity can still invalidate this approximation.
    """
    selected = found[(found['size']>=low)&(found['size']<=high)]
    blocks = (selected.start_bin.to_numpy()*dt//block_ms).astype(int)
    total_blocks = int(np.ceil(span/block_ms))
    hist = np.zeros((total_blocks, high-low+1),dtype=np.int64)
    np.add.at(hist,(blocks,selected['size'].to_numpy().astype(int)-low),1)
    counts = hist.sum(axis=0)
    fit,lp = fit_counts(counts,low)
    residual = np.cumsum(counts/counts.sum())-np.cumsum(np.exp(lp))
    # All-zero blocks can be aggregated without altering the multinomial law.
    nonzero = np.any(hist,axis=1)
    nonempty = hist[nonzero]
    probabilities = np.r_[np.full(len(nonempty),1/total_blocks),1-len(nonempty)/total_blocks]
    weights = rng.multinomial(total_blocks,probabilities,size=reps)[:,:-1]
    boot_counts = weights.astype(float) @ nonempty.astype(float)
    support = np.arange(low,high+1)
    taus, boot_lp = scalar_family(boot_counts,np.log(support))
    boot_residual = np.cumsum(boot_counts/boot_counts.sum(axis=1)[:,None],axis=1)-np.cumsum(np.exp(boot_lp),axis=1)
    centered = np.max(np.abs(boot_residual-residual),axis=1)
    return dict(block_ms=block_ms,block_tau_low=float(np.quantile(taus,.025)),
                block_tau_high=float(np.quantile(taus,.975)),
                block_ks_p=float((1+np.count_nonzero(centered>=fit['ks']))/(reps+1)),
                block_count=total_blocks)
