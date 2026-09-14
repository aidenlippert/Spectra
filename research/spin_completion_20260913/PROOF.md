# What spin completion proves

All claims concern the frozen rational H6 electronic Hamiltonian, with twelve
spin orbitals and six electrons. The ten retained spatial density matrices
and the exact collective remainder are unchanged. The retained Hamiltonian
H10 obeys H10 + ell I <= H, where ell = -3.632e-9 Ha. No new molecular
integrals or many-body wavefunction enter the lower discovery.

## Additional joint operators

For the same rational matrices L_k, define

    Q_k^(s,t) = sum_pq L_kpq a^dagger_(p,s) a_(q,t).

The previous frame used only s=t. The completed frame admits both diagonal
and spin-changing components, and all products Q_k^(s,t) a_i together with
a_i. Individual spin-changing components need not be Hermitian. The exact
CAR polynomial is therefore formed first; its adjoint is computed from that
polynomial, reversing operator order correctly.

The 492 listed generators partition into four 93-dimensional and four
30-dimensional blocks by exact spin charge and molecular parity. Each old
252-frame generator is literally included. The numerical map tests check
this containment and both new spin-charge and mixed-component moment maps.
No numerical basis truncation is performed.

Whitening is only an invertible numerical congruence used for conditioning.
The additional full-family control omits it, preserving every listed
generator and every PSD entry. Both coordinate choices export into the same
original exact polynomial basis. The conditioning choice therefore changes
numerical work and possible solver error, not the mathematical cone.

For any selected real rational combinations B_j, a PSD matrix G supplies
sum_jk G_jk {B_j^dagger, B_k} >= 0. This includes cross-combination terms.
Restricting G to a nonnegative diagonal gives the separate-contribution
control in exactly the same B basis.

## Exact lower acceptance

The compact certificate stores two integer coordinate matrices: selected
directions V over the frame, and a rounded Gram factor R over those selected
directions. Compose R V as integers with the product denominator. Expand
each resulting operator D exactly, then construct D^dagger exactly. The
ordinary CAR verifier receives both squares D^dagger D and D D^dagger.
Their degree-six terms cancel; acceptance refuses any surviving term above
degree four. This catches numerical solve, map, rounding, and export errors.

With the unchanged full quadratic baseline and Hermitian body-one X,

    H10 - b I = sum T^dagger T + sum {D^dagger,D}
               + (Nhat-6)X + Rres.

The exact number ideal vanishes on the fixed-N sector and each fermionic
monomial has norm at most one. Therefore b - ||Rres||_1 + ell is a valid
lower for the original H. The separate rational Rayleigh witness provides
an upper for that same H. A numerical solver status is not an acceptance
criterion. Only exact reconstruction and replay certify each interval.

## Exact counterexample separation versus energy gain

The preceding full diagonal-spin dual y is independently replayed against
all its original PSD, number-ideal, coefficient-box and binding conditions.
A separator is an integer frame vector v for which the exactly evaluated
y({B^dagger,B}) is strictly negative. This proves that the additional
positivity condition excludes that particular nonphysical functional.

A negative direction alone proves no bound gain and no new optimum. Energy
improvement comes from separately optimized, exactly accepted lower proofs.
If a new lower exceeds the previous complete-frame dual ceiling, it proves
that the improvement cannot be reproduced anywhere in the previous cone.
Failure to reach 1.6 mHa does not by itself establish a new obstruction.

## Cost and scope

The first adaptive pricing round uses the accepted counterexample's exact
separators; later rounds use the current numerical dual only as a proposal.
Both tested bundle sizes, shared construction, exact acceptance, rejected
candidates and exports count toward the outer adaptive watchdog. Prior
counterexample discovery and its source solve are recorded as ancestry.
The diagonal control inherits all selected-direction discovery cost.

The ten Hamiltonian patterns and the selected correlation combinations are
different objects. Neither count is a count of many-body configurations or
an asymptotic compression theorem. The stored upper's previous determinant
amplitude discovery remains outside the compact lower method. This finite
basis test does not establish general many-body solvability or physical
preparation and control of matter.
