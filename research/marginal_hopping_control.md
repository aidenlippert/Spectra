# Full mixed-dictionary hopping control

I ran the requested six-mode, (N=3) controls at (t=1/5), with one matched and one asymmetric inter-triple hopping pattern. Each solve used

    baseline_blocks() + seed_blocks() + cubic_pools("mixed")

This is an expressivity ceiling for the finite sector, not a scalable speed baseline. The proposal solves took 11.53 s (matched) and 6.30 s (asymmetric), with one thread and the configured Python/CVXPY stack.

| case | certified lower | variational upper | certified width | residual row norm |
|---|---:|---:|---:|---:|
| matched | 0.550995396532 | 0.551000400320 | 5.00e-6 | 5.246082e-6 |
| asymmetric | 0.540863615947 | 0.540877181134 | 1.36e-5 | 1.363602e-5 |

Receipts: [matched](../results/marginal_hopping/full_mixed_matched.json), [asymmetric](../results/marginal_hopping/full_mixed_asymmetric.json), [summary](../results/marginal_hopping/full_mixed_summary.json).

## Verifier audit

The verifier reconstructs each factor as an integer linear combination of exact occupation-basis word matrices, forms (\sum F^\dagger F\), and compares the residual to the exact rational Hamiltonian. It checks Hermiticity before estimating the residual. The reported `residual_row_norm` is the induced infinity norm (maximum absolute row sum). For a symmetric/Hermitian residual, this is a valid upper bound on the spectral norm, so the returned lower bound (b-\eta) is safe. The residual check is stronger than trusting the SDP status or objective.

The trial vector is an integer occupation-basis vector and the upper bound is computed as the exact rational Rayleigh quotient (x^\dagger Hx/(x^\dagger x)). Thus the upper number is variational and does not depend on the floating-point reference eigenvalue; `reference_numeric` is diagnostic only.

The export still clips tiny negative floating Gram eigenvalues before factorization and rounds factors to a fixed decimal denominator. This can enlarge the residual, which is correctly charged to (\eta). It does not establish a production-scale rational PSD decomposition, and the full mixed dictionary remains a finite-sector ceiling. The result should therefore be read as: exact replay certifies a (5\times10^{-6}) to (1.4\times10^{-5}) energy bracket for these two instances, not as evidence of transfer or scaling.

