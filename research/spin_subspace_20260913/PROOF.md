# Selected spin combinations and exact energy bounds

This pass keeps the same fixed-N=6 H6 Hamiltonian and ten retained spatial
density patterns. All energies are electronic Hartree in the frozen rational
finite-basis model. The [previous collective proof](/Users/aidenlippert/Documents/Spectra/research/molecular_collective_20260913/PROOF.md)
supplies the exact operator statement

\[
H_{10}+\ell I\preceq H\preceq H_{10}+uI,
\quad \ell=-3.632\times10^{-9},\quad u=1.08327632\times10^{-4}.
\]

## Two small coordinate maps, then exact adjunction

Write Q_k=Q_k,up+Q_k,down using the same rational spatial matrix L_k for both
spins. In each exact molecular symmetry class form a frame Phi of a_i and
Q_k,spin a_i. There are four frames of dimension 63. These are operators,
not many-electron configurations.

The learned combinations are B_j=sum_i V_ji Phi_i, where the coordinate
matrix V is rationalized before solving. Their joint positive contribution
has the form sum_jk G_jk {B_j^dagger,B_k}, with G positive semidefinite.
A real numerical Gram factor is rounded to a rational matrix R. Each proof
row is therefore a linear combination of the original frame with coefficient
row (R V)_r.

The compact certificate stores the integer numerators of R and V and their
two positive denominators. Acceptance composes the maps as integer dot
products, expands each resulting row D_r exactly, and constructs its partner
D_r^dagger by exact adjunction. Its contribution is

\[
D_r^\dagger D_r+D_r D_r^\dagger\succeq0.
\]

No independently rounded adjoint is admitted. The existing CAR SOS verifier
then expands both squares. The degree-six leading terms cancel exactly;
acceptance refuses a surviving degree-six residual. Numerical map errors,
PSD clipping, direction rounding, factor rounding, and solve errors are all
exposed by this exact reconstruction.

The baseline remains the full quadratic fermionic cone. With Hermitian
one-body X, a certificate gives

\[
H_{10}-bI=\sum T^\dagger T+\sum_r\{D_r^\dagger,D_r\}
+(\widehat N-6)X+R_{\mathrm{res}}.
\]

The number ideal vanishes in the chosen sector and each CAR monomial has
norm at most one. Hence the accepted molecular lower is
b-||R_res||_coefficient,1+ell. The rational reference upper is separately
replayed against H itself, so the interval uses the original Hamiltonian
at both endpoints.

## Why the ablation measures cross-combination contributions

The joint solve permits arbitrary PSD G in the selected B basis. The separate
solve permits only diagonal nonnegative G in that same basis. Its feasible
cone is contained in the joint cone. Export may collect its separate square
rows in a common block without adding cross terms: the Gram matrix is still
the sum of those row outer products. Every direction and its normalization
is identical between the two controls.

The measured difference between the two exact lower certificates establishes
a gain for the found joint proof. It does not, by itself, give an exact
optimum or obstruction for the separate cone.

## What direction selection establishes

The full-frame moment matrices are computed from the current numerical dual.
Their negative eigenvectors are only proposals. Numerical independence is
tested in the operator coefficient space, and directions are rounded before
their candidate optimization. Every accepted gain comes from exact lower
replay. The four- and eight-direction bundles are compared using their actual
gain divided by their measured construction/solve/export/replay cost plus
the round's pricing cost.

All candidate solves, including those not selected for continuation, remain
charged in the total adaptive budget. Selecting by a measured rate does not
establish an advantage over another ranking policy without a matched policy
comparison. In particular, choosing eight directions in every round supplies
no evidence of superiority to a fixed eight-direction schedule.

## A full-frame ceiling applies to every learned subspace

Let y be a Hermitian linear functional with y(I)=1, |y(w)|<=1, and zero
expectation for every one-body number-ideal polynomial. Suppose its quadratic
moment matrices and all four complete spin-frame anticommutator matrices
are positive semidefinite. Then all selected moment matrices are positive:

\[
\Gamma_B=V\Gamma_\Phi V^T\succeq0.
\]

Consequently y(H10)+ell is a ceiling for every lower certificate made from
any such subspace, including any future choice or number of directions inside
these same frames. This follows directly by applying y to the certificate
identity and bounding y(R_res)>=-||R_res||_1. It needs no assumption of
numerical convergence or strong duality.

The exact full-frame checker first reuses the previous binding, coefficient,
normalization, Hermiticity, symmetry, number-ideal, and quadratic PSD gates.
It then checks all four 63-dimensional matrices explicitly. The old
spin-summed patterns are sums of the two spin components and therefore lie
inside these frames. Their additional inherited PSD checks are redundant,
not a substitution for the full spin checks. The test suite verifies that
the exact checker rebuilds every generator used by numerical discovery and
rejects the previously feasible spin-summed witness.

An accepted full-frame dual ceiling C can be compared with the independently
replayed prior full-cubic physical lower L*: for every present lower L,
E0(H)-L>=L*-C. This bounds physical error only if that comparison is made
against the same frozen Hamiltonian and particle sector. The audit verifies
the prior polynomial, sector, certificate hash, and spectral residual proof.

All of these statements concern the specified quadratic-plus-frame cone,
one-body number ideal, and coefficient-L1 residual rule. They do not exclude
other correlated operators, spin-changing density components, different
residual certificates, or other compact many-body representations.
