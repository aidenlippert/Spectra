# Hamiltonian-only interaction commutant route

The useful structural question is whether the quartic interaction itself
reveals a hidden density basis. Let (V) be the number-conserving quartic
part and let (Q(A)=\sum_{pq}A_{pq}a_p^\dagger a_q) for a Hermitian one-body
matrix (A). Form the exact linear map

\[
\mathcal C_V:A\longmapsto[V,Q(A)].
\]

This is a linear system with (M^2) unknowns and, after normal ordering, at
most (O(M^4)) word equations. Its kernel is computable without a supplied
orbital basis or source factors.

## Recovery theorem

Suppose the interaction is density-density in an unknown one-particle basis:

\[
V=\sum_{i,j} U_{ij} n_i n_j+\sum_i u_i n_i,
\qquad n_i=d_i^\dagger d_i,
\]

and suppose the kernel of (mathcal C_V) contains an (M)-dimensional real
self-adjoint pairwise-commuting linear space (mathcal A). Any self-adjoint
basis of (mathcal A) can be simultaneously diagonalized. Because there are
M linearly independent commuting matrices on an M-dimensional one-particle
space, require that their joint diagonal map has rank M; it then spans the
full diagonal space and contains the rank-one orbital projectors as linear
combinations. Their joint eigenvectors recover the (d_i) modes up to phases
and permutation. Transforming (V) by
that recovered unitary gives a density-density coefficient tensor, which can
then be tested for additional path, locality, or stoquastic structure.

The proof is elementary: pairwise commuting self-adjoint matrices are simultaneously
unitarily diagonalizable, and rank M of their diagonal image gives every
rank-one projector. Since the interaction commutes with every
recovered number operator, its normal-ordered quartic terms preserve each
occupation label and therefore have density-density form in that basis. The
kernel need not be closed under multiplication; this argument does not assume
it is an algebra.

The hypotheses matter. Kernel dimension (M) alone is insufficient: the
kernel may be noncommutative, its commuting diagonal image may have rank below
M, or it may contain
accidental symmetries unrelated to number operators. Degenerate (U) can
produce a larger nonabelian commutant, in which case the basis is not
identifiable without (H_1), higher-body terms, or a gauge convention.

## Exact and numerical implementation boundary

For rational CAR coefficients, construct (mathcal C_V) over exact rationals,
compute its nullspace, and explicitly verify pairwise commutators and
self-adjointness before diagonalization. Also require a rank-one-projector or
joint-spectrum gate; repeated joint eigenvalues do not identify orbitals.
Numerical eigendecomposition is suitable for
discovery only; an irrational or approximate rotation cannot be inserted into
an exact certificate without an outward error bound. A rational orthogonal
rotation, when present, can be replayed exactly. Otherwise certify a nearby
orthogonal U0 by interval eigenspaces or rational orthogonal approximation.
If both U and U0 are exactly orthogonal/unitary and ||U-U0||op <= delta, a
telescoping bound gives the certified transformation error
4 delta ||V_4||1 + 2 delta ||V_2||1, where V4 and V2 are the quartic and
one-body coefficient L1 masses. If the proposed matrix is not exactly unitary,
it cannot transport the CAR algebra; first replace it by an exact orthogonal
matrix and certify that replacement error.

The proposed workflow is therefore:

1. extract (V) from (H), build the exact commutator map, and compute its
   kernel;
2. find an M-dimensional pairwise-commuting self-adjoint linear space, verify
   that its diagonal image has rank M, and recover the rank-one projectors;
3. recover the orbital rotation, with exact or certified interval data;
4. transform (H), test density/path structure, and only then emit SOS factors;
5. replay the transformed certificate and bound the rotation/truncation error.

Exact nullspace computation has polynomial equation count but bit complexity
depends on coefficient heights; numerical recovery additionally needs a
nonzero joint eigengap. This route is strongest for deliberately rotated
density-density families.
Generic molecular interactions need not have an (M)-dimensional commutant;
spin (SU(2)), spatial, or point-group symmetries commonly give only a small
nonabelian symmetry algebra. Thus commutant recovery is a Hamiltonian-only
structure detector, not a general molecular solver or a consequence of
low-rank tensor representations.
