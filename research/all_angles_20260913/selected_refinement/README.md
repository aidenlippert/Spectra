# Signed-residual selected CI: a new upper-witness path for the finite ladder

This is an implementation of established selected configuration interaction principles, not a new general many-body algorithm. It was built during the multi-agent campaign after the first dense CISD/selected-CI probes missed H8 accuracy.

The generator reads the original rational Hamiltonian and a supplied single-determinant Hartree–Fock reference. It compiles each CAR word into occupation/parity masks with integer coefficients and combines cancellations exactly for each requested determinant action. A selected subspace grows by ranking the signed residual coupling to omitted determinants, divided by a regularized diagonal-energy difference. SciPy's sparse eigensolver optimizes the state within the selected subspace. It does not enumerate the full Fock space up front, read saved FCI coefficients, or read a lower certificate during candidate generation.

Each selected eigenvector is rounded to integer amplitudes. The existing standard-library streaming checker computes its exact Rayleigh quotient for the original H. Only that quotient is an upper endpoint. A small projected eigensolver residual establishes neither ground-state identity nor a many-body error bar.

`replay_ladder.py` checks original-H and sector equality, hashes the old lower certificate and new upper witness, then replays the complete lower proof and the new upper state. Existing lower certificates and, for H6/H8/H10, existing spectral residual witnesses are reused. Their historical discovery cost is not removed by this work.

The first passing selected supports are 16, 128, 2,048, and 24,576 for H4, H6, H8, and H10 respectively. This rapid growth is evidence against claiming scalable accuracy from this ladder. H4/H6 were run locally; H8/H10 were run on an existing Lambda A10 host using one CPU BLAS thread. GPU hardware was not used by this sparse algorithm. The separate GPU track tested a different projected eigensolve.

The search parameters were fixed support ladders, an amplitude/candidate proposal threshold of 1e-14, a denominator floor of 1e-3, and integer amplitude rounding at 1e12. Those approximations can weaken the trial state but cannot invalidate the exact upper replay. All intermediate witnesses, including intervals too wide to pass, are retained.

Timing in each discovery row is cumulative through that row and includes its exact upper replay. It excludes the historical lower-certificate search and the separate complete interval replay. The complete H10 initial 16,384-state run and the subsequent longer rerun are separate charges; the extended run did not secretly inherit its cache. Final report accounting includes both. Peak memory across the whole run was not instrumented; sparse entries, action-cache size, candidate counts, and logical term/state checks are recorded instead.

Focused tests compare compiled action against the independent CAR action on several H4 determinants, replay rounded states exactly, and check invalid-sector and budget refusal. The independent audit also inspected the mathematical distinction between a numerical proposal and a certified upper endpoint.

Example, from the repository root with NumPy and SciPy available:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.all_angles_20260913.selected_refinement.refine \
  --fixture results/certificate_scaling/active_space_ladder/h8/fixture.json \
  --reference results/certificate_scaling/active_space_ladder/h8/upper.json \
  --out results/all_angles_20260913/selected_refinement/new_h8_run \
  --budgets 32 64 128 256 512 1024 2048 4096 --seconds 180
```

The generator refuses an existing output directory. Certificate replay uses `python -S` via `replay_ladder.py`; no floating-point library is needed for checking the frozen endpoints. Exact arithmetic does not prove the Python implementation formally correct or address continuum/basis/model errors.
