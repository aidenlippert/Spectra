# SCS resume design for adaptive PSD discovery

The current `adaptive_block_discovery.run` rebuilds a CVXPY `Problem` and all
variables on every refinement round, then calls `problem.solve(solver="SCS",
...)` without `warm_start=True`. Its `save_raw` files contain the exported
free variables, residuals, dual rows, and PSD values, but not SCS's canonical
cone state (`x`, `y`, `s`). They cannot resume an SCS iterate directly.

CVXPY's installed SCS interface does support warm starts. In the local
CVXPY 1.9.2 source (`SCS.solve_via_data`), when `warm_start` is true and a
`solver_cache` entry exists, it passes cached `x`, `y`, and `s` to `scs.solve`.
After a status of exactly `optimal`, it stores the returned result in that
cache. Thus the smallest safe seam for repeated solves of an unchanged
canonical problem is to retain one `Problem` object and call:

```python
problem.solve(solver="SCS", warm_start=True, **options)
```

The cache is process-local. A resumable artifact would need to save the
canonical solver inputs and state explicitly: sparse `A`, vectors `b` and `c`,
cone dimensions, SCS version, options, and `x/y/s`. Store arrays in NPZ and a
JSON manifest containing SHA-256 hashes, matrix dimensions/nnz, cone data,
solver/formulation/conditioning parameters, and the exact support list. On
load, reconstruct the canonical data and refuse resume unless every hash,
dimension, cone field, option, and SCS version matches. This prevents applying
an iterate to a changed support set or differently scaled problem.

For adaptive refinement, adding supports changes the cone dimensions and
canonical `A`; an old `s` is then not valid. The safe seam is therefore:

1. Keep the same `Problem` and solver cache for bounded repeated solves before
   pricing, using `warm_start=True`.
2. When supports change, rebuild and hash a new canonical problem. Reuse only
   mapped primal/dual values if an explicit dimension-preserving embedding is
   implemented and checked; otherwise start the new cone solve cold.
3. Save the new canonical state only after SCS returns a usable result, with
   `optimal` distinguished from `optimal_inaccurate`.

The remote warm hosts use CVXPY 1.6.5 and Clarabel 0.11.1 in
`/home/ubuntu/spectra-venv`; local is CVXPY 1.9.2. The CVXPY SCS cache behavior
should be verified against the remote installed source before relying on
cross-version serialized state. No evidence currently shows that
`optimal_inaccurate` results are cached by CVXPY: the 1.9.2 source caches only
when status equals the exact `optimal` token. A direct `scs.solve` wrapper is
needed if bounded continuation from an inaccurate iterate is required.

The existing H10 observation (175 SCS iterations taking roughly 673 seconds)
therefore motivates keeping one process and warm-starting identical repeated
solves, but it does not establish a safe cross-round resume until canonical
data identity and cone-shape checks are implemented. No timing or speedup is
claimed here.

## Implemented bounded canonical restart

Both installed interfaces were inspected directly: local CVXPY1.9.2/SCS3.2.11
and remote CVXPY1.6.5/SCS3.2.8. Both cache only optimal results, not
optimal_inaccurate. Merely requesting warm_start=True cannot repair that.
The earlier statement about an absent warm_start argument must not be read
as a claim that CVXPY's default is false; the decisive issue is the cache.

scs_checkpoint.py now calls the conic solver interface with an explicit
canonical x/y/s cache, including finite status-2 iterates. It hashes A/P/b/c,
cone data, versions, relevant solver settings, and the original problem
contract. Changed time/iteration budgets are allowed; changed mathematics,
versions or other options are refused. Array files are hashed and loaded
without pickle; shapes and finite float64 values are checked. Existing
checkpoint files are preserved. It does not restore SCS's internal adaptive
scaling, acceleration history, or factorization; this is an iterate restart,
not an identical uninterrupted trajectory.

Two real-SDP tests cover inaccurate status persistence, restart into a fresh
CVXPY Problem, convergence, changed budgets, changed data/contracts, and byte
corruption refusal. They passed locally and on both remote environments.
The full H4 pipeline resumed after two deliberately short 50-iteration stages
and exported an exact lower certificate; its separate final interval is in
results/certificate_scaling/scs_refinement/h4_resumed/interval.json.

Two H10 jobs compare one1200-second solve with two600-second restarted solves.
They have equal nominal solver time, while setup/restart/export costs remain
separately charged. Their statuses are live experimental state, not results;
inspect the named handles and remote processes before reporting completion.
