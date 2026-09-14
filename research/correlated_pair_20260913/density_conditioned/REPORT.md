> Superseded exploratory report. Do not use its conclusions as accepted evidence. See [the final campaign report](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/REPORT.md) and [rejected prototypes](/Users/aidenlippert/Documents/Spectra/research/correlated_pair_20260913/REJECTED_PROTOTYPES.md).

# Density-conditioned pair branch

This branch implements two concrete pieces of the correlated-factor mechanism.

The closure builder keeps

    F = psi_i + n_j sum_k d_k psi_k + tau + mu n_j
        + epsilon n_j nu,

with cubic `nu` commuting with `n_j`. Exact symbolic multiplication retains
the quartic terms and density-times-cubic terms of degree six. The checker
rejects a closure that drops those degree-six terms. A four-mode resonant
example produces such terms, so quartic-only closure is not accepted.

The second component is an exact rational max-flow feasibility checker for the
restricted parent family H0 = sum_i e_i n_i + sum_ij U_ij n_i n_j. For every
negative edge it seeks endpoint allocations `w_ij,i + w_ij,j = -U_ij`, with
node load at most `e_i`. Positive U terms are separate positive squares. The
max-flow/min-cut result is an exact certificate or infeasibility witness.

The H6/H8 diagnostic first performs the particle-hole transform around the
Hartree-Fock occupation pattern (modes `0..N-1`), expands the diagonal density
polynomial exactly, and adds the smallest rational chemical-potential shift
that makes every transformed linear capacity positive. It then runs the exact
flow/min-cut test on the transformed pair coefficients. This avoids the trivial
failure from feeding signed bare molecular one-body coefficients directly into
the parent rule. The resulting feasibility or min-cut witness is recorded in
`results/.../diagnostic.json`; it remains a restricted-family result, not a
molecular-spectrum certificate.

All arithmetic in the checker is rational; no FCI or Gram solve is used.
