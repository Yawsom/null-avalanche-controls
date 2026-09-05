# Reproduce the feasibility audit

Read [REPORT.md](REPORT.md) for the decision, [PROTOCOL.md](PROTOCOL.md) for the
experimental design, and [LITERATURE.md](LITERATURE.md) for source comparisons.
The next proposed study is described in [RESEARCH_PLAN.md](RESEARCH_PLAN.md);
its experiments have not yet been run.

Use Python 3.12. The recorded environment used Python 3.12.14 on Windows.
The pinned packages in `requirements-lock.txt` include the original repo's
dependencies and SciPy for independent/constrained optimization. The original
`requirements.txt` remains unchanged.

From the repository root, in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r feasibility/requirements-lock.txt
$env:MPLCONFIGDIR = Join-Path $PWD '.mplcache'
$env:OPENBLAS_NUM_THREADS = '1'
.\.venv\Scripts\python.exe bp_validate.py
.\.venv\Scripts\python.exe -m unittest feasibility.test_methods -v
.\.venv\Scripts\python.exe -m feasibility.run --stage all
.\.venv\Scripts\python.exe -m feasibility.summarize
.\.venv\Scripts\python.exe -m feasibility.verify
```

On Linux/macOS, use `python3 -m venv .venv`, the executable `.venv/bin/python`,
and `export MPLCONFIGDIR="$PWD/.mplcache" OPENBLAS_NUM_THREADS=1`.
The Windows platform was tested; the other platforms were not.

Stages `main`, `long`, `controls` and `calibration` can also run independently.
They write separate result files and print progress. A completed condition is
checkpointed to CSV, but a stage restart recomputes that stage from the start.
Allow several minutes and additional time for installation; the main grid took
about 4.5 minutes in the audited environment. Wall time depends on hardware and
concurrent workloads. No large raw event files are saved.

`summarize` requires all four stages to be complete. It regenerates the summary
tables and PNG/SVG figures, including the seed-42 distribution plot. The report
is a written interpretation of the recorded run; it is not automatically
rewritten if experiment parameters change.

`verify` checks source hashes and result completeness against this protocol. It
also checks the archived validation logs. To refresh those logs, redirect the
two validation commands' stdout/stderr to `feasibility/results/validation_legacy.txt`
and `feasibility/results/validation_feasibility.txt` respectively. Source-hash
checks are expected to fail if the numerical code is edited after a run.
Verification accepts LF/CRLF conversion and either manifest path separator so
Git checkout conventions do not invalidate otherwise identical source files.

The manifests capture simulation/statistical source hashes, seeds, dependencies
and the baseline Git commit. They record the state when each stage began;
nonexecuted documentation/test/plotting files can be added subsequently. The
original numerical source files and the two experiment kernels (`run.py` and
`methods.py`) should match the manifests used for their results.

The original scripts are intentionally retained as the baseline. The new
detector wrapper validates inputs; the new finite-array generator is a
sensitivity model, not an in-place replacement advertised as critical. The
statistical routines do not claim that KS acceptance or surrogate rejection
establishes criticality or causation.
