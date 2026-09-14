# H6 cubic certificate precision diagnosis

Audited read-only artifact:
`results/lambda_runs/structural_scaling/downloaded/cubic_small/h6_mixed_full/`.
The certificate has 12 modes, 6 particles, 58 retained supports, 1,108
factor rows, 121,994 factor nonzeros, 9,419 square-polynomial terms, and a
maximum residual degree of six. The solver history reports one round,
`optimal_inaccurate`, numeric lower `-6.333617717744763`, exact lower
`-6.3410288521244`, and dual identity `0.9999999999999807`.

Re-expanding the saved rational certificate with the exact CAR verifier gives
the following residual decomposition (absolute coefficient l1, term count,
largest coefficient):

| body degree | l1 | terms | max abs coefficient |
|---:|---:|---:|---:|
| 0 | 0.00275960253104 | 1 | 0.00275960253104 |
| 2 | 0.00057629411791 | 36 | 0.00004618252088 |
| 4 | 0.00322996810838 | 882 | 0.00002089309841 |
| 6 | 0.00408873636707 | 8,500 | 0.00000668082839 |

These sum to the recorded exact residual l1 `0.0106546011244`. The saved
multiplier has 919 terms, coefficient l1 approximately `6.137715302`, and
maximum absolute coefficient approximately `1.039569754`. The factor matrices
are available only after integer factor export; their aggregate count is
121,994 nonzeros. The run cost was 7.70 s map construction, 61.26 s solve,
10.45 s exact export, 79.43 s total.

The dominant residual contribution is degree six (38.4% of l1), followed by
degree four (30.3%), degree zero (25.9%), and degree two (5.4%). This shows
that the loss is distributed across the generated CAR hierarchy rather than
being visibly isolated to a single low-degree equality.

The rounded certificate alone cannot identify whether the error came from PSD
eigenvalue clipping, solver feasibility error, factor quantization, or
rationalization of `b` and `X`. The exported integer factors determine a
valid exact Gram, but they do not preserve the pre-export floating Gram
eigenvalues or the clipped directions. The `optimal_inaccurate` status and the
gap between numeric and exact lower are evidence of numerical sensitivity, not
a causal attribution. Raw PSD blocks/eigenvalues, pre-rounding residuals, and
solver primal/dual residuals are required to separate those mechanisms.
