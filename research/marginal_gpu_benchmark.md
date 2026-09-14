# Frozen restricted LP backend benchmark

`results/marginal_graded_hubbard8/discovery/gpu_lp_benchmark.py` rebuilds a
frozen numerical restricted LP from the accepted joint Reynolds crossover
proposal. It preserves proposal ordering, writes the matrix through HiGHS as
an MPS file, and records row residuals and objective value from a CPU HiGHS
run. Passing `--cuopt` attempts the same MPS input only when an importable
cuOpt installation is already present; no package installation or cloud
provisioning is performed.

Run from the repository root:

```text
python -m results.marginal_graded_hubbard8.discovery.gpu_lp_benchmark \
  --seconds 120 --out results/marginal_graded_hubbard8/gpu_benchmark
```

Add `--cuopt` for an optional cuOpt probe. The receipt distinguishes an
unavailable package, API failure, and a completed backend result. This is a
numerical comparison of a frozen restricted model. It does not create or
validate a polynomial certificate; exact rational export and independent
replay remain the acceptance gate.

The corrected CPU fixture contains 1,239 rows, 6,660 columns (4,148 atom directions plus signed metric coefficients and residual variables), and 874,793 nonzeros. Native interior point plus crossover reached objective approximately −9.09e−14 with maximum equality residual 2.52e−11 in 27.83 seconds. An earlier invalid fixture omitted metric variables and produced the trivial objective 1.001; it is not evidence for this workload.

The optional GPU branch now uses the documented `Problem.readMPS`, `SolverSettings`, `Variable.Value` and `Constraint.DualValue` interfaces. It restores primal/dual ordering by MPS names and records residuals, bound violations and solution arrays. `--cuopt-method` selects barrier or PDLP. This branch has not run on a GPU. See the [NVIDIA API reference](https://docs.nvidia.com/cuopt/user-guide/latest/cuopt-python/lp-qp-milp/lp-qp-milp-api.html).

Lambda credentials were not found in the current environment, standard shell/project configuration, or standard Lambda configuration locations. A credential-location question is pending. No instance was created and no credits were spent.
