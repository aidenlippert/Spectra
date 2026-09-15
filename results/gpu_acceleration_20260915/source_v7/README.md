# GPU acceleration experiment

This directory is isolated from the main transfer campaign. It contains the
measured CPU/CUDA proposer and the unchanged-checker replay wrapper. No main
solver or accepting-checker source was edited.

The selected backend sends blocks of size 128 or larger to an A100 FP64
eigensolver and keeps smaller blocks on the CPU. Ordinary batched CUDA eigh
was slower on these fixtures. `cpu_evd` is the CPU fallback.

The separate `--adaptive` option changes the numerical search schedule. Starting
at mu=2, it switches to mu=.03 after the predicted interval falls below 10 mHa
and the primal residual falls below 1e-4. Both are checked every 100 iterations.
The existing numerical stopping requirements (predicted interval below 1.3 mHa
and primal residual below 5e-8) remain. Neither prediction accepts a certificate:
the original rational checker must establish the complete interval below
1.6 mHa. This schedule change is not a GPU-only speedup.

## Run a prepared case

From the repository root, using an environment with NumPy, SciPy, threadpoolctl,
and CuPy plus its CUDA toolkit dependencies:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m research.gpu_acceleration_20260915.solve CASE RESULT --seconds 240 --mu 2 --backend hybrid --adaptive
python -B -S -m research.gpu_acceleration_20260915.replay CASE RESULT
```

`CASE` must contain the preparation, fixture, integer MPS, upper receipt, and
nonsinglet proof expected by the transfer solver. `RESULT` must be new.
The proposer validates the fixture binding and normal-map identity. Replay
checks both the MPS upper and the complete-spin lower with exact arithmetic;
it refuses changed input hashes or inconsistent endpoints. GPU proposals are
never accepted using floating-point eigenvalues alone.

`--iterations N` provides an iteration-matched benchmark instead of early
termination. Use the same schedule, inputs, precision, and host for a hardware
comparison. `kernel_bench.py` additionally measures synchronized projection
timings with and without transfer costs.

## Fresh construction

`cold_run.py NEW_CASE --state-python /absolute/path/to/state_venv/bin/python`
constructs the asymmetric H6 fixture from new integrals and runs every step.
It records source hashes before starting, gives every stage a time limit, and
preserves failures. It inherits no earlier state, coefficient maps, optimization
checkpoint, or proof. Software installation and machine provisioning are
separate costs; CUDA library caches may already be warm.

The tested state environment uses Python 3.12, quimb 1.15.0, NumPy 2.2.6,
SciPy 1.18.1, networkx and opt_einsum. The CUDA numerical environment uses
Python 3.10, NumPy 2.2.6, SciPy 1.15.3, CuPy 14.2.0, and threadpoolctl 3.6.0.
The complete construction additionally requires PySCF, numba, CVXPY, and
SCS 3.2.11. Exact installed versions are retained with the run receipts.

## Verification and exploratory branches

```sh
OPENBLAS_NUM_THREADS=1 python -B -m unittest research.gpu_acceleration_20260915.test_projection
OPENBLAS_NUM_THREADS=1 SPECTRA_TEST_GPU=1 python -B -m unittest research.gpu_acceleration_20260915.test_projection
```

Projection tests include indefinite, degenerate, zero and negative matrices,
symmetry, positivity, CPU/GPU agreement, and invalid input. Complete solver
outputs are also checked against CPU checkpoints and independently replayed.

`exact_arithmetic_probe.py` is an experimental GMP comparison, not the retained
accepting path. It requires an existing authoritative stdlib receipt and
demands exactly equal rational endpoints. Its timings varied, and GMP integer
transfers did not demonstrate an improvement. It is not enabled by default.

See `REPORT.md` for measured results and the failed environment/setup probes.
