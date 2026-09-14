# Fresh H10 lower-bound audit (read-only feasibility)

The required original-H10 fixture is
`results/certificate_scaling/active_space_ladder/h10/fixture.json`, SHA256
`ec437b52d210f08233a8cc853d0f905dae52bd1f2a923285a25c76ac99385f92`.

The old baseline is `results/lambda_runs/scs_refinement/downloaded/h10_evd_1800`.
Its command (from `h10_evd_1800.log`) was:

```text
/home/ubuntu/spectra-venv/bin/python research/certificate_scaling/adaptive_block_discovery.py --fixture results/certificate_scaling/active_space_ladder/h10/fixture.json --outputdir results/h10_evd_1800 --family mixed --ideal-body 2 --symmetry --full --solver SCS --solver-eps 1e-8 --solver-seconds 1800 --seconds 2300 --solver-max-iters 150000 --denominator 1000000000 --free-denominator 1000000000000 --save-raw --quotient-linear --anchor-identity --condition row-column --save-scs-state
```

It is a fresh SOS discovery from H, with no source factors or upper witness.
Measured map construction was 154.227 s, SCS solve 1817.291 s, exact export
251.174 s, total 2222.948 s. The exact coefficient-L1 lower was
`-12.376685805916855 Ha`, residual L1 `0.007797540375855352 Ha`, with 4,674,900
complete pricing entries, 2,347,238 map nonzeros, 4,810 factor rows, and
2,581,696 factor nonzeros. The archived certificate hash is recorded in the
campaign report; this is an historical baseline, not a new endpoint.

The fresh run target is lower `>= -12.3698994564 Ha`, corresponding to the
existing fresh upper minus 0.0016 Ha. A bounded reproduction should use the
same command and exact fixture hash, with output isolated under
`results/wave2_20260913/fresh_lower/h10_fresh_*`; do not substitute
`spin_irrep/*/symmetric_hamiltonian.json`.

Dependencies: `/home/ubuntu/spectra-venv/bin/python`, CVXPY/SCS, repository
`research/certificate_scaling` sources, and (for the old speed path) the EVD
shim/library if enabled by the host. Expected resources are roughly 4.4 GB
resident memory and 37 minutes for the old 1800-second protocol. No execution
was performed in this read-only audit pending Lambda A assignment.

## Smaller fresh spin-invariant baseline

Run on the exact original-H10 fixture:

```text
/home/ubuntu/spectra-venv/bin/python research/certificate_scaling/spin_invariant_discovery.py --fixture results/certificate_scaling/active_space_ladder/h10/fixture.json --out results/wave2_20260913/fresh_lower/h10_spin_scs600 --solver SCS --seconds 600
```

Then generate and independently replay the spectral witness:

```text
/home/ubuntu/spectra-venv/bin/python research/certificate_scaling/wedge_spectral_bound.py --certificate results/wave2_20260913/fresh_lower/h10_spin_scs600/certificate.json --out results/wave2_20260913/fresh_lower/h10_spin_scs600_spectral
```

Both receipts are required; solver status alone is not an endpoint. This is
fresh discovery with original H and no source factors. Root's matched H6 run
showed only a `2.693877e-6 Ha` gain from full-number-X, so this is a baseline
measurement rather than an assumed cure.
