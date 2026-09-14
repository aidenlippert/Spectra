# A directly generated local positivity proof

The Hubbard-chain complement certificate can now be constructed and replayed without a configuration list, numerical optimizer, or full-population number-ideal completion. Its expanded positive atoms have occupation degree at most four and number O(m²) for m sites. Exact replay passes through 32 spatial sites, the existing 64-mode verifier limit.

This is a deliberately loose complement lower bound for a structured Hamiltonian. It establishes a complete bounded-degree construction for this family. It does not establish accurate energy intervals, molecular transfer, or efficient general marginal representability.

## Formula and exact decomposition

Take an even m, fixed Nα=Nβ=m/2, and a nearest-neighbor chain with L=m−1 edges:

\[
H=U\sum_i n_{i\alpha}n_{i\beta}
-t\sum_{\langle i,j\rangle,\sigma}
(a^\dagger_{i\sigma}a_{j\sigma}+a^\dagger_{j\sigma}a_{i\sigma}),
\quad U,t\ge0.
\]

Let D count doublons and Q be the ionic complement D≥1. Choose any rational ε>0 and

\[
v=\tfrac12\sum_i q_i^2,
\qquad \gamma=U-4tL-\varepsilon.
\]

On the fixed-number slice, v=D. Globally on Boolean occupations, v=D−T/2, with T=Nα+Nβ−m. Hence the metric proof is the exact polynomial identity

\[
v=1+(D-1)-\tfrac12(N_\alpha-m/2)-\tfrac12(N_\beta-m/2).
\]

For one edge and one spin, write If=n_jσ(1−n_iσ), Ib=n_iσ(1−n_jσ). A forward hop changes D by n_i,barσ−n_j,barσ. Define

\[
S=(1-I_f-I_b)(D+1)
+I_f(1-n_{i\bar\sigma}+n_{j\bar\sigma})
+I_b(1-n_{j\bar\sigma}+n_{i\bar\sigma}).
\]

This is nonnegative on all Boolean assignments: 1−If−Ib equals n_iσ n_jσ+(1−n_iσ)(1−n_jσ), and each remaining factor is a sum of nonnegative occupation indicators. Expanding S requires only indicators of degree at most four. Contradictory occupation events vanish; duplicate atoms combine.

Summing over the 2L edge/spin pairs gives the exact identity

\[
K_D=(UD-\gamma)D-t\sum_{\text{directed hops}}I D'
=U D(D-1)+(2tL+\varepsilon)(D-1)+\varepsilon+t\sum S.
\]

Every displayed nonconstant contribution is nonnegative on Q. The actual compiler uses v=D−T/2, so its numerator differs by the explicit number identity

\[
K_v=K_D-\tfrac{T}{2}
\left(UD-\gamma-t\sum_{\text{directed hops}}I\right).
\]

Both spin-number multipliers are therefore the degree-two polynomial −(UD−γ−tΣI)/2. No lifting to full-population monomials is needed. Existing coefficient replay verifies the metric lower bound 1 and numerator lower bound ε with exactly zero residual. It reconstructs the encoded Hamiltonian and target, rather than trusting a formula label.

Since v(0)=0, transitions to P carry zero weight. The projector-free weighted-row theorem then certifies QHQ≥γ. The bound is not an accuracy improvement over elementary row estimates; for example, unweighted row counting gives the stronger coarse chain bound U−2tL. The purpose here is a complete local-atom construction without the sector-sized machinery required by the optimized H6 proof.

## Measured size transfer

The following runs use U=4, t=1/3, ε=1/1000. Determinant actions, sector generation, and the full-population ideal-lift helper are patched to reject calls during construction and measured replay. Each saved certificate is then replayed independently with Python's standard library only.

| Sites | Spin-sector dimension, not enumerated | Positive atoms | Number-multiplier terms | Certificate bytes | Exact replay seconds |
|---:|---:|---:|---:|---:|---:|
| 4 | 36 | 72 | 38 | 13,776 | 0.0082 |
| 6 | 400 | 160 | 58 | 25,615 | 0.0178 |
| 8 | 4,900 | 280 | 78 | 40,796 | 0.0279 |
| 10 | 63,504 | 432 | 98 | 59,380 | 0.0451 |
| 16 | 165,636,900 | 1,080 | 158 | 137,389 | 0.1839 |
| 24 | 7,312,459,672,336 | 2,392 | 238 | 298,509 | 0.6760 |
| 32 | 361,297,635,242,552,100 | 4,216 | 318 | 533,264 | 1.9417 |

At 32 sites the constructor takes approximately 0.021 seconds; replay performs 100,921 polynomial coefficient products. Timings are individual measurements on the current machine, not a general runtime guarantee. The atom count is quadratic, number-identity terms are linear, and coefficient degree stays bounded. Generic exact spin-symmetry checks and other runtime costs are included in replay timing.

The 32-site complement bound is −112003/3000. There is no corresponding ground-energy interval in this result. The mathematical formula holds for all even m; the implementation deliberately retains the existing 32-site cap, denominator budgets, exact-sign checks and other verifier refusal paths.

## Validation and remaining target

Three focused tests exercise exact transfer at 2, 4, 6, 8 and 10 sites; all physical CAR rows at four sites including P sources; U=0 and t=0; invalid parameters; and corruptions of the Hamiltonian, threshold and positivity coefficients. The accepted proof path is `experiments.marginal_polynomial_metric.replay`, shared with the molecular proofs. Every acceptance requires coefficient reconstruction and exact positivity checks against the encoded Hamiltonian.

Artifacts: `results/marginal_hubbard_polynomial/{sites}/certificate.json`, `receipt.json`, `independent_replay.json`, and `scaling.json`. Constructor: `experiments/marginal_hubbard_polynomial.py`.

The unresolved step is to combine comparably controlled construction cost with the accuracy of the optimized H6 certificates and a compressed retained-response calculation. [H6 results and remaining costs](marginal_joint_metric.md).

Full marginal validation: **353 tests pass**; seven new chain certificates independently replay with the standard library.
