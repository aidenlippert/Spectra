# Measured optional libraries for Spectra

This isolated package calls the existing proposal solver with sparse QR and
optional cuDSS/CHOLMOD preconditioning. All experiments are separate processes;
the frozen solver and accepting-source files are preserved.

Use a prepared case containing `fixture.json`, `upper.json`, `nonsinglet.json`
and the existing `prepared` coefficient arrays. The recorded Hamiltonian and
normal-map consistency checks still apply. A new result name is required.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m research.nvidia_followup_20260915.solve CASE RESULT --seconds 850 --backend hybrid --normal cudss
python -B -S -m research.nvidia_followup_20260915.strict_replay CASE RESULT NEW_EXACT_DIRECTORY
```

For the H10 orbital-rotation upper, append `--rotated ROTATED_H10_DIRECTORY` to
replay. The public replay caller enforces the paired-mode, even-electron and
nonempty-nonsinglet assumptions of this spin argument. An odd-electron problem
requires a different spin-sector construction. This explicit guard does not
alter any accepted campaign fixture, which all satisfy those assumptions.

The CPU proposal configuration is `--backend cpu_evd --normal superlu`.
`--normal cholmod` selects the CPU sparse Cholesky alternative. The optional
FLINT rational replay uses `--compiled-rationals` and an environment containing
python-flint (omit `-S` in that case). Its retained experimental results also
require an independent standard-library replay with matching exact endpoints.

The remote environment uses Python 3.12, NumPy 2.2.6, SciPy 1.18.1, CuPy 14.2,
nvmath-python 1.0, python-flint 0.9, sparseqr 1.6 and system SuiteSparse 5.10.1.
The full installed-version list is in the downloaded environment record.
`bootstrap.sh` and `accelerated.sh` record installation commands for the isolated
Lambda environment; they are not intended to change a local global environment.

Focused checks:

```sh
python -B -m unittest research.nvidia_followup_20260915.test_sparse_quotient research.nvidia_followup_20260915.test_flint_squares
python -B -S -m unittest research.nvidia_followup_20260915.test_sector_gate
```

`library_bench.py` measures actual sparse normal matrices and an actual MPS
transfer. `accelerated_jobs.py` records iteration-matched solver comparisons
before a fresh integrated H10 optimization. `replay.py` supplies the existing
upper and lower checks behind the public sector guard. Float predictions are
never accepted as certificates.

The additional exact integer-Gram implementation was validated but was not
selected, because its measured improvement was small. cuTensorNet was likewise
not selected for the tested small MPS transfer. See REPORT.md for measured
effects, scope, complete dependency costs and the terminated cloud-instance
receipt. LIBRARIES.md contains the broader chemistry-software shortlist.
