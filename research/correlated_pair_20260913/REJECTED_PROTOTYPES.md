# Rejected exploratory implementations

The early `self_consistent/builder.py`, `density_conditioned/maxflow.py`,
`density_conditioned/closure.py`, and `mps_search/mps_dmrg.py` are retained
as development history and deliberately refuse execution. Their preliminary
reports and tests do not support the final campaign claims.

Integration found an incorrect creator-map prefactor and extra weight,
incorrect treatment of the Hermitian conjugate half and scalar term,
particle-hole density/sign errors, a same-sign chemical shift, incomplete
flow witnesses, a concatenation surrogate for CAR products, and an incomplete
MPS expectation/export adapter. These were implementation errors, not
mathematical obstructions. No final result uses those implementations.

The accepted replacements are `self_consistent/fixed_guide.py`,
`density_conditioned/exact_density.py`, `mps_direct.py`, `mps_spatial.py`,
`mps_round.py`, and `mps_exact.py`. The root campaign report supersedes the
early branch reports. Final validation uses `test_mps`, `test_algebra`,
`self_consistent.test_fixed_guide`, and `test_tensor_adapter` only.

Two early failed optimizer output files are empty and have no usable timing
receipt. Their logs remain in `results/correlated_pair_20260913/runs`.
The final cost account explicitly marks these times unavailable.
