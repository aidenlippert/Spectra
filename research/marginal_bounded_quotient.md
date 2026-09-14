# Degree-preserving coefficient quotients and positivity transfer

The accurate H6 certificate now has a proof-construction path that avoids both physical-configuration evaluation and lifting monomials to full particle population. It uses exact fixed-number coefficient quotients, then checks the exported certificate with the existing Hamiltonian-bound rational verifier. The metric and some positive directions are inherited from earlier finite-sector discovery; that historical cost has not disappeared.

## The algebraic replacement

For m Boolean variables and fixed population p, let V_d contain multilinear occupation polynomials of degree at most d. Assume d≤min(p,m−p). The relations

\[
g_S=(\sum_i x_i-p)x_S,\qquad |S|\le d-1,
\]

are reduced by x_i²=x_i. They span exactly the polynomials in V_d that vanish on the p-particle slice. The quotient has dimension C(m,d), rather than sum_{k≤d} C(m,k). This dimension is polynomial in m for fixed d; it says nothing about which degree is sufficient for accurate positivity certificates.

The implementation eliminates these relations using exact rational arithmetic and descending (degree,mask) order. Every elimination row carries its explicit signed number multiplier. Consequently the output is an exact coefficient identity f=r+(N−p)L, with degree(r)≤degree(f) and degree(L)≤degree(f)−1. Alpha and beta reductions are applied separately while carrying the opposite-spin factors. No occupation assignments are generated.

The previous helper raised monomials until each spin population was fully occupied. The new helper stops within the input degree. On H6 the degree already reaches the three-particle boundary of each spin sector, so the joint quotient still has 20×20=400 coordinates and the reconstructed proof still has 1,364 number-multiplier terms. A different coordinate system alone is not compression of this particular full function space.

## Why the quotient is complete in the stated degree range

On homogeneous subset coefficients, let U add one element to a subset and D=U* remove one. On degree k,

\[
DU-UD=(m-2k)I,
\qquad \|Uf\|^2=\|Df\|^2+(m-2k)\|f\|^2.
\]

Thus U is injective for k<m/2. The highest-degree part of (N−p)L is U applied to the highest-degree part of L. Since degree(L)≤d−1<m/2, the displayed number relations are linearly independent.

For p≤m/2, evaluating degree-d monomials on p-element subsets is the injective map U^(p−d)/(p−d)!. Lower-degree monomials can be raised only to degree d using

\[
x_S=\frac{\sum_{i\notin S}x_{S\cup\{i\}}}{p-|S|}
\]

on the slice. The evaluation rank of V_d is therefore C(m,d). Subtracting this from dim(V_d) gives exactly the number of independent relations, proving they span the evaluation kernel. For p>m/2, the substitution x_i→1−x_i preserves degree and reduces evaluation rank to the complementary slice. The same relation-independence argument applies.

The code checks the expected exact rank as well as its degree and monomial budgets. The existing positivity verifier independently reconstructs every exported identity; quotient completeness is not a new trusted acceptance shortcut.

## Fresh coefficient-space transfer

For the existing 1/50 hopping perturbation, the numerical constructor uses the original-H vanishing metric, discards the original positive-atom weights, and keeps its positive directions as a seed. It adds every positive indicator of degree at most four that has more than one fixed-spin completion, plus every nonzero degree-four-or-lower charge localizer. Completion counts are combinatorial formulas, not configuration enumeration.

The resulting LP has 400 quotient-coordinate rows, 20,110 columns and 1,480,387 nonzeros. Preparation took approximately 22.86 seconds and the solve 15.94 seconds. It returned a positive numerator bound of approximately 0.01370188943. Five tiny negative proposed coefficients, the largest in magnitude 1.74e−8, were clipped explicitly during proposal export. The final rational residual gate accepted the reconstructed certificate; LP status and tolerance did not supply the acceptance.

The fresh perturbed numerator uses 299 positive indicators and 93 charge localizers. The metric lower bound after the new number-identity reconstruction is 996951691/10^12. Both the original and perturbed proof reconstructions retain all 1,364 number-multiplier terms. They preserve the old metric and require no new full-population lifting.

This is new positivity-coefficient discovery for the changed Hamiltonian, with an inherited metric and inherited seed directions. It is not a Hamiltonian-only discovery algorithm, a proof of the entire marginal cone, or a claim that all 400 H6 coefficient coordinates have been compressed away.

## Degree-capped scaling diagnostic

Exact coefficient identities were reconstructed without evaluating physical assignments:

| One-spin variables | Degree cap | Input monomials | Number relations | Quotient dimension | Physical-slice dimension, not enumerated |
|---:|---:|---:|---:|---:|---:|
| 6 | 2 | 22 | 7 | 15 | 20 |
| 8 | 3 | 93 | 37 | 56 | 70 |
| 12 | 3 | 299 | 79 | 220 | 924 |
| 16 | 3 | 697 | 137 | 560 | 12,870 |
| 24 | 2 | 301 | 25 | 276 | 2,704,156 |
| 32 | 2 | 529 | 33 | 496 | 601,080,390 |

The 32-variable degree-two quotient constructs in approximately 0.0037 seconds in this run. These are coefficient-algebra measurements, not large-system positivity or energy certificates. The implementation preserves a 4,096-monomial construction budget and refuses work outside its supported degree range.

## Evidence and remaining work

Source: `experiments/marginal_number_quotient.py`; focused tests: `tests/test_marginal_number_quotient.py`. The tests reconstruct exact coefficient identities, check every assignment of several small symmetric and asymmetric slices, exercise the two-spin identity, and test invalid inputs and degree/budget refusals.

Artifacts and reproducible scripts live under `results/marginal_h6/polynomial_metric/bounded_quotient/` and its parent `discovery/` directory. The machine-readable progress record distinguishes independent replays, numerical proposals, inherited work and unresolved outcomes.

The remaining targets are controlled metric/direction discovery on larger molecular systems, useful accuracy at bounded positivity degree, and a retained-response representation that avoids explicit many-body configurations. [Earlier joint-metric results](marginal_joint_metric.md) and [the analytic local family](marginal_hubbard_polynomial.md) retain their separate scopes.

Final validation: **360 full marginal tests pass**. The current progress record accounts for all nine new independent standard-library replays.
