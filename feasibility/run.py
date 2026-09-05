"""Run with python -m feasibility.run --stage all (see PROTOCOL.md)."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import scipy
from scipy.special import gammaln

from bp_analyze import average_iei, DEFAULT_BINS_MS
from .methods import (absolute_fit, block_fit, generate, interval_test, measure,
                      metrics, size_counts, surrogate)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'feasibility'/'results'
CONDITIONS = ['homogeneous','uniform','mixed','capped','refractory_0.8',
              'refractory_1.0','refractory_1.2']
SEEDS = list(range(100,120))


def save(rows, filename):
    pd.DataFrame(rows).to_csv(OUT/filename,index=False)


def signature(row, estimator):
    return (abs(row['sigma_single']-1)<.2 and abs(row[estimator]-1.5)<.2
            and row.get('vuong_z',np.nan)>2)


def run_record(condition,hours,seed,bootstrap=199):
    _,times,meta = generate(condition,hours,seed)
    base = dict(condition=condition,hours=hours,seed=seed)
    iei = average_iei(times,200)
    native = max(1,int(round(iei)))
    widths = sorted(set(map(int,DEFAULT_BINS_MS.split(',')))|{native})
    sweep, cache = [],{}
    for dt in widths:
        found = measure(times,dt)
        cache[dt] = found
        row = dict(**base,dt_ms=dt,iei_ms=iei,native_dt_ms=native,
                   **metrics(found))
        row['signature_regression'] = signature(row,'tau_regression')
        row['signature_mle'] = signature(row,'tau_mle') if 'tau_mle' in row else False
        sweep.append(row)
    frame = pd.DataFrame(sweep)
    best_reg = float(frame.loc[np.hypot(frame.sigma_single-1,frame.tau_regression-1.5).idxmin(),'dt_ms'])
    best_mle = float(frame.loc[np.hypot(frame.sigma_single-1,frame.tau_mle-1.5).idxmin(),'dt_ms'])
    policies = {'iei_rounded':float(native),'fixed_10':10.,
                'selected_regression':best_reg,'selected_mle':best_mle}
    detail, memo = [],{}
    for policy,dt in policies.items():
        found = cache[dt]
        for low,high in ([(1,30),(1,15),(1,59),(3,30)] if policy in ['iei_rounded','fixed_10'] else [(1,30)]):
            key=(dt,low,high)
            if key not in memo:
                m=metrics(found,low,high)
                counts,_=size_counts(found['size'].to_numpy(),low,high)
                if counts.sum()>=50:
                    rng=np.random.default_rng(np.random.SeedSequence([seed,int(hours*100),int(dt),low,high,777]))
                    m.update(absolute_fit(counts,low,rng,bootstrap))
                memo[key]=m
            row=dict(**base,policy=policy,dt_ms=dt,iei_ms=iei,
                     fit_min=low,fit_max=high,**memo[key])
            row['signature_regression']=signature(row,'tau_regression')
            row['signature_mle']=signature(row,'tau_mle') if 'tau_mle' in row else False
            detail.append(row)
    # Separate numerical/bin-rule sensitivity; do not treat it as a new tuned rule.
    sensitivity=[]
    for label,dt in [('iei_float',iei),('iei_100',max(1,round(average_iei(times,100)))),
                     ('iei_400',max(1,round(average_iei(times,400))))]:
        sensitivity.append(dict(**base,policy=label,dt_ms=dt,**metrics(measure(times,dt))))
    if condition in ['capped','refractory_0.8','refractory_1.0','refractory_1.2']:
        f=cache[4]
        meta.update(detected_fraction_over_20ms=float(np.mean(f.duration*4>20)),
                    detected_max_duration_ms=int(f.duration.max()*4))
    return sweep,detail,sensitivity,dict(**base,iei_ms=iei,**meta)


def main_grid():
    all_sweep,all_detail,all_sensitivity,all_meta=[],[],[],[]
    start=time.monotonic()
    for hours in [.25,3.]:
        for condition in CONDITIONS:
            for seed in SEEDS:
                s,d,x,m=run_record(condition,hours,seed)
                all_sweep.extend(s);all_detail.extend(d);all_sensitivity.extend(x);all_meta.append(m)
            save(all_sweep,'sweeps.csv');save(all_detail,'fits.csv')
            save(all_sensitivity,'iei_sensitivity.csv');save(all_meta,'generators.csv')
            print(f'main hours={hours} condition={condition} complete; {time.monotonic()-start:.1f}s',flush=True)


def long_grid():
    sweeps,details,sens,meta=[],[],[],[]
    for seed in [42,100,101,102,103,104]:
        s,d,x,m=run_record('mixed',70.,seed)
        sweeps.extend(s);details.extend(d);sens.extend(x);meta.append(m)
        save(sweeps,'long_sweeps.csv');save(details,'long_fits.csv')
        save(sens,'long_iei_sensitivity.csv');save(meta,'long_generators.csv')
        print(f'long seed={seed} complete',flush=True)
    blocks=[]
    # Dependence-aware sensitivity at 3h and 70h; no pooled event-level intervals.
    for condition,hours,seed in [('mixed',3.,100),('capped',3.,100),
                                  ('refractory_1.0',3.,100),('mixed',70.,42)]:
        _,times,meta=generate(condition,hours,seed)
        native=max(1,round(average_iei(times,200)))
        for dt in sorted({native,10,4}):
            found=measure(times,dt)
            for block_ms in ([1000,10000] if hours==3. else [10000]):
                result=block_fit(found,dt,meta['span_ms'],np.random.default_rng(seed+dt+block_ms),block_ms=block_ms)
                blocks.append(dict(condition=condition,hours=hours,seed=seed,dt_ms=dt,**result))
        save(blocks,'block_bootstrap.csv')
        print(f'block sensitivity {condition} {hours}h complete',flush=True)


def controls():
    paired,intervals=[],[]
    for ci,condition in enumerate(CONDITIONS):
        for seed in SEEDS:
            _,times,meta=generate(condition,.25,seed)
            iei=average_iei(times,200); native=max(1,round(iei))
            for kind,width in [('none',0),('shuffle',0),('jitter',4),('jitter',20),('jitter',80)]:
                rng=np.random.default_rng(np.random.SeedSequence([seed,ci,width,7777]))
                moved=times if kind=='none' else surrogate(times,rng,meta['span_ms'],kind,width)
                new_iei=average_iei(moved,200); new_native=max(1,round(new_iei))
                for dt in sorted({native,10,new_native}):
                    paired.append(dict(condition=condition,seed=seed,control=kind,jitter_ms=width,
                                       dt_ms=dt,original_native=native,control_native=new_native,
                                       **metrics(measure(moved,dt))))
            for width in [20,40,100]:
                rng=np.random.default_rng(np.random.SeedSequence([seed,ci,width,8888]))
                intervals.append(dict(condition=condition,seed=seed,window_ms=width,
                                      **interval_test(times,rng,meta['span_ms'],width)))
        save(paired,'surrogate_metrics.csv');save(intervals,'interval_tests.csv')
        print(f'controls {condition} complete',flush=True)


def calibration():
    rows=[]
    s=np.arange(1,31,dtype=float)
    pmfs={'exact_powerlaw':s**-1.5,
          'borel_critical':np.exp(-s+(s-1)*np.log(s)-gammaln(s+1)),
          'exponential':np.exp(-.2*s)}
    for label,p in pmfs.items():
        p=p/p.sum()
        for n in [1000,10000]:
            for seed in range(200,250):
                rng=np.random.default_rng(seed)
                counts=rng.multinomial(n,p)
                fit=absolute_fit(counts,1,rng)
                rows.append(dict(distribution=label,n=n,seed=seed,**fit))
        print(f'calibration {label} complete',flush=True)
    save(rows,'calibration.csv')


def manifest(stage):
    paths=[ROOT/x for x in ['bp_simulate.py','bp_analyze.py','bp_validate.py','simulate.py']]
    paths+=list((ROOT/'feasibility').glob('*.py'))+[ROOT/'feasibility'/'PROTOCOL.md']
    record=dict(stage=stage,python=sys.version,platform=platform.platform(),
                numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,
                base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seeds=SEEDS,hours=[.25,3.],bootstrap_replicates=199,interval_replicates=99,
                conditions=CONDITIONS)
    (OUT/f'manifest_{stage}.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['all','main','long','controls','calibration'],default='all')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    stages={'main':main_grid,'long':long_grid,'controls':controls,'calibration':calibration}
    for stage,action in stages.items():
        if args.stage in ['all',stage]:
            manifest(stage);action()


if __name__=='__main__':
    main()
