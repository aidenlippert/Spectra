# Conditioning polynomial dictionaries without changing their operators

The H6 Krylov failures motivated a numerical experiment on the same exact
polynomial span. This is not a certificate pruning theorem.

Write the original generators as a column p, and their monomial coefficient
matrix as C (monomials by generators). A square basis change p'=W*p gives
C'=C*W-transpose. A PSD Gram Q in the new coordinates represents the old
Gram W-transpose*Q*W. The contraction backend uses C' throughout its numerical
map. Export factors Q=Z-transpose*Z are multiplied as Z*W, then rounded into
integer coefficients of the ORIGINAL rational polynomials. Thus every exported
factor stays in the original span, even when the numerical transformation has
rounding error. The production verifier pays the complete physical residual.
No numerical invertibility assertion is used as an energy certificate.

Whitening computes a thin SVD and retains every direction. Numerically rank-
deficient bases are refused rather than automatically truncated. The implementation
first identifies connected components of generators sharing exact monomial
support, then whitens each component separately. This keeps off-component
zeros exactly zero. The orthogonality check and SVD timing are recorded.

Disjoint coefficient supports do NOT make physical operators independent.
Cross-component Gram entries remain in the SDP and tests explicitly check this.
The component graph is used only to avoid needless numerical fill in W.

## H4 measurements

On the identical diagonal depth-one, creator-channel H4 dictionary:

| Basis | Map nonzeros | Build s | Solve wrapper s | Total s | Exact interval Ha |
|---|---:|---:|---:|---:|---:|
| Original contraction | 167104 | .681 | 3.395 | 4.466 | .0000651119 |
| Whole-block SVD | 2392701 | .618 | 20.309 | 21.183 | .00000121250 |
| Exact-support component SVD | 167104 | .604 | 2.393 | 3.235 | .00000123319 |

All intervals have independent python -S replays. Whole-block SVD introduced
floating fill even across exactly orthogonal support components. The component
version avoids it without threshold clipping. The whole-block experiment is
retained as evidence of this performance pitfall, not deleted as an abandoned
run. These are single measurements, not a robust asymptotic speedup estimate.

## H6 measurements and scope

The component transformation keeps all 43740 Gram entries and 9419 coefficient
rows. Its proposal map has 4839646 nonzeros (two fewer floating nonzeros than
the previous map; no coefficient threshold is applied). Construction is 11.64
seconds. In the first linear block, the exact-support components have sizes
1,1,1,6,6,6,6,6,4,6,6,6,4,6,4; its numerical orthogonality error is 1.83e-14.

The Clarabel120 experiment finished at its time limit. It took 137.53 seconds
in the solve wrapper and 161.71 seconds overall. Its exact lower is
-6.389140221365018 and its interval is approximately .0560816 Ha, still above
.0016. This improves the .353744 interval of the earlier unwhitened Clarabel
run but is not a passing chemistry result or a proof of the cone optimum.
The matched SCS180 run also completed and independently replayed: width
.210586928014 Ha, construction11.53s, solve231.22s, total254.29s. It is
also a failure against .0016, despite improvement over4.24923 unwhitened.

Raw numerical objective, negative eigenvalue mass, raw b, exact coefficient
residual, export time and factor sizes are retained in each receipt. Negative
eigenvalue mass is coordinate-dependent and must not be compared as though
it were an invariant physical error. Whitening does not remove the polynomial
products, coefficient rows, certificate bit growth or physical reference cost.

An audit suggestion that a diagonal commutator can never add a direction was
incorrect. It assigns a different energy weight to each monomial, generally
changing the polynomial's direction. A new exact test gives a two-monomial
counterexample. H4 and H6 numerical basis ranks were full at the tested
threshold; H6 was not shown to have more near-dependence than H4. Likewise,
for equation scaling D, physical dual values are D*y_scaled, as implemented.

Changed modules: polynomial_gram_contraction.py, commutator_dictionary.py.
37 focused tests pass, including congruence orientation, every transformed
CAR-map coefficient on a nontrivial control, rank refusal, retained cross-
component variables and an end-to-end positive-Hamiltonian solve. No full
repository regression or generic structural accuracy theorem is claimed.

A depth-zero H4 control, using the same conditioning, remains above the
known exact obstruction: independently replayed width .003566801899 Ha.
Whitening has not bypassed that restricted cone's mathematical limitation.
