"""Verify result completeness and the numerical source hashes used in each stage."""
import hashlib
import json
from pathlib import Path
import pandas as pd


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'feasibility'/'results'
    checks={}
    names={'run.py','methods.py','bp_simulate.py','bp_analyze.py','bp_validate.py','simulate.py'}
    for stage in ['main','long','controls','calibration']:
        manifest=json.loads((out/f'manifest_{stage}.json').read_text())
        kernels={k.replace('\\','/'):v for k,v in manifest['source_sha256'].items()
                 if k.replace('\\','/').split('/')[-1] in names}
        def matches_recorded_hash(path, expected):
            raw=path.read_bytes()
            lf=raw.replace(b'\r\n',b'\n')
            variants=(raw,lf,lf.replace(b'\n',b'\r\n'))
            return any(hashlib.sha256(value).hexdigest()==expected for value in variants)
        checks[stage+'_kernel_hashes']=all(matches_recorded_hash(root/k,v)
                                           for k,v in kernels.items())
    sweep=pd.read_csv(out/'sweeps.csv');fits=pd.read_csv(out/'fits.csv')
    long=pd.read_csv(out/'long_fits.csv');interval=pd.read_csv(out/'interval_tests.csv')
    recordings=sweep[['condition','hours','seed']].drop_duplicates()
    checks['280_main_recordings']=len(recordings)==280
    checks['20_seeds_each_cell']=bool((recordings.groupby(['condition','hours']).size()==20).all())
    checks['no_duplicate_sweeps']=not sweep.duplicated(['condition','hours','seed','dt_ms']).any()
    checks['2800_fit_rows']=len(fits)==2800
    checks['60_long_fit_rows']=len(long)==60
    checks['420_interval_tests']=len(interval)==420
    checks['300_calibration_samples']=len(pd.read_csv(out/'calibration.csv'))==300
    checks['no_lognormal_bounds']=not fits.lognormal_bound_hit.any() and not long.lognormal_bound_hit.any()
    checks['bootstrap_p_range']=bool(fits.ks_iid_p.dropna().between(.005,1).all())
    checks['legacy_checks_pass']='all checks passed' in (out/'validation_legacy.txt').read_text()
    log=(out/'validation_feasibility.txt').read_text()
    checks['eight_new_tests_pass']='Ran 8 tests' in log and 'OK' in log
    (out/'verification.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(checks,indent=2))
    if not all(checks.values()):raise SystemExit('Verification failed')


if __name__=='__main__':main()
