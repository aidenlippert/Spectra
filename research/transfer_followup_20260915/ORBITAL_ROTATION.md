# Exact orbital rotation for a compact upper certificate

This construction changes the representation of a trial state while retaining
an upper bound for the original supplied rational Hamiltonian. It does not
assume that a numerical orbital transformation or an approximate MPO is exact.

Let the spatial transformation be \(V=Z/d\), where \(Z\) is an integer square
matrix, \(d>0\), and the checker verifies

\[
Z^T Z=d^2 I.
\]

Apply the same transformation to both spins. The spin-orbital transformation
\(O=V\otimes I_2\) is exactly orthogonal. Its second quantization is a unitary
\(\Gamma(O)\) on Fock space that preserves particle number and total spin.
The checker never constructs this full many-body matrix.

Write the normal-ordered Hamiltonian as

\[
H=cI+\sum_{pq}h_{pq}a_p^\dagger a_q
  +\sum_{p<q,r<s}C_{pq,rs}a_p^\dagger a_q^\dagger a_s a_r.
\]

Its exact transform has coefficients

\[
h'=O^T hO,\qquad
C'=(\wedge^2 O)^T C(\wedge^2 O),
\]

where
\((\wedge^2 O)_{pq,ij}=O_{pi}O_{qj}-O_{pj}O_{qi}\).
These are one- and two-orbital operator coordinates. For 20 spin orbitals the
two-orbital array has 190 coordinates; it is not the fixed-10-electron space.
The implementation further partitions these coordinates by alpha count, which
the paired spatial rotation preserves.

The proposal uses a numerically designed local orbital basis. Rational
Householder reflections approximate that design while preserving orthogonality
exactly: for an integer vector \(v\),

\[
R_v=I-2vv^T/(v^Tv),\qquad R_v^TR_v=I.
\]

The approximation to the desired local basis affects efficiency, not validity:
any exactly orthogonal accepted rotation gives the unitary identity above.

Round the transformed polynomial to a rational grid to obtain \(\widetilde H\).
The checker computes every rounding difference exactly and sums

\[
\epsilon=\sum_w |h'_w-\widetilde h_w|.
\]

Every fermionic monomial has operator norm at most one. Consequently,
\(\|\Gamma(O)^\dagger H\Gamma(O)-\widetilde H\|\le\epsilon\).
For any nonzero, exactly fixed-number integer MPS \(\psi\), the original MPS
checker evaluates its rational Rayleigh quotient \(U_{\rm rot}\). Therefore

\[
E_0(H\vert_N)\le
\frac{\langle\psi,\Gamma(O)^\dagger H\Gamma(O)\psi\rangle}
{\langle\psi,\psi\rangle}
\le U_{\rm rot}+\epsilon.
\]

This is the accepted original-model upper. Numerical DMRG sweep energies,
approximate MPO contraction energies, and DMRG convergence flags are not used
as accepted endpoints. A completed-sweep checkpoint can be checked even when
the numerical optimizer has not declared convergence.

The paired lower calculation keeps the original canonical Hamiltonian and its
existing full fixed-number SOS checker. A newly constructed HF product supplies
only proposal moments. It need not be the state supplying the accepted upper.
The final interval combines that lower with the separately checked rotated-state
upper and charges the discovery and verification of both components.

The focused test compares the transformed one- and two-body polynomial against
an independently implemented bit-state matrix on a four-mode example, verifies
the rounding inequality, and rejects a nonorthogonal transformation. This is
an algebraic implementation check, not an independent full H10 oracle.

There is no claim that all orbital transformations, MPS bond dimensions, or
induced lower representations remain cheap at arbitrary size. This experiment
tests a specific measured state-construction bottleneck.
