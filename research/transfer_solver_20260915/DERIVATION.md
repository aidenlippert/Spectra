# Input-derived spin certificates and the scope of the transfer test

This pass generalizes the successful H8 construction to dimensions and spatial
parities derived from the supplied Hamiltonian. It retains the earlier exact
accepting implementation. The numerical constructor and the mathematical
certificate have different responsibilities: an optimizer proposes integer
factors; exact operator arithmetic decides their energy guarantee.

## Hamiltonian and target

The input is an explicitly supplied rational, Hermitian, number-conserving
operator on paired alpha/beta spin orbitals,

\[
 H=\sum_w h_w w,\qquad h_w\in\mathbb Q,\qquad
 \{a_p,a_q^\dagger\}=\delta_{pq},\quad
 \{a_p,a_q\}=0.
\]

The accepted interval concerns the full fixed-even-particle-number sector
\(\mathcal H_N\). The molecular fixtures contain electronic energies; adding
the same nuclear-repulsion constant to both endpoints leaves the width unchanged.
The target is a **total** interval width of at most \(1/625\) hartree, or
1.6 mHa. It is not a per-electron target.

Fresh floating-point molecular integrals are rounded to an explicit rational
Hamiltonian. The recorded coefficient norm encloses this rounding relative to
those floating inputs. It does not enclose integral evaluation, basis error,
geometry error, or a difference between the mathematical model and experiment.

## Exact spin reduction and the two lower bounds

Let \(\mathcal T(A)=\int_{SU(2)}U_g^\dagger A U_g\,dg\) denote normalized
spin averaging, and write

\[
 H_s=\mathcal T(H),\qquad
 \delta=\sum_w|[H-H_s]_w|.
\]

Every creation or annihilation operator has norm at most one. Consequently,
\(\|H-H_s\|\le\delta\). The original checker constructs the average in
exact rational CAR arithmetic and checks that both spin ladder commutators
of \(H_s\) vanish. It does not assume that rounded integrals have exact spin
symmetry.

The singlet proof has the form

\[
 H_s-bI=\mathcal T\!\left(\sum_j A_j^\dagger A_j+J+R\right),
\]

where every \(A_j\) is an integer-coefficient operator polynomial divided by
a declared integer denominator. The accepted ideal term satisfies
\(P_0 J P_0=0\), with \(P_0\) the singlet projector within \(\mathcal H_N\).
Its allowed components are

\[
 (\widehat N-N)X,
 \quad(\widehat N_\alpha-N/2)Y,
 \quad a\widehat S^2,
 \quad S_+W+(S_+W)^\dagger.
\]

The checker enforces the degree, Hermiticity and conservation requirements on
these multipliers. In particular, \(W\) lowers \(M_S\) by one. Each component
has zero singlet compression. Spin averaging preserves positivity of each
square, so

\[
 L_0=b-\sum_w|[\mathcal T(R)]_w|
 \le \lambda_{\min}(H_s|_{S=0}).
\]

A separately constructed certificate proves a bound \(L_1\) on the complete
\(M_S=1\) sector, without using singlet-only identities there. For even
particle number, every total-spin multiplet has integer spin. Every
\(S\ge1\) multiplet has an \(M_S=1\) component with the same energy for
the spin-invariant operator \(H_s\). Therefore

\[
 \boxed{L=\min(L_0,L_1)-\delta\le\lambda_{\min}(H|_{\mathcal H_N}).}
\]

The defect \(\delta\) is charged once. Neither an \(M_S=0\) trial nor an
\(M_S=0\) numerical eigenvector is treated as a proof of pure singlet spin.
Odd-electron models and empty nonsinglet screens require a different sector
rule; this constructor refuses those cases.

## What is compressed by the highest-weight representation

The mixed operator dictionary contains linear and cubic fermionic words.
It is partitioned by particle-number change, twice spin projection and
input-derived binary spatial parity. Only parities commuting with spin
rotation are used. Both spin-ladder closure conditions are checked, and the
raising matrices are compared with literal CAR commutators in focused tests.

Within each parity and charge class, the mixed dictionary contains spin
\(1/2\) and spin \(3/2\) representations, including all their multiplicities.
In normalized irreducible-tensor coordinates, Haar orthogonality gives

\[
 \mathcal T(O_{j,m,\alpha}^\dagger O_{j',m',\beta})
 =\frac{\delta_{jj'}\delta_{mm'}}{2j+1}
   \sum_{q=-j}^{j}O_{j,q,\alpha}^\dagger O_{j,q,\beta}.
\]

Thus the average of an arbitrary positive Gram matrix depends on positive
matrices in the **multiplicity** coordinates. The trace over magnetic
coordinates remains positive. Conversely, any such positive multiplicity
matrix can be represented by positive squares built from one highest-weight
component. This is an equivalence of the full spin-averaged cones, provided
all multiplicities are retained. An invertible change to rational,
unnormalized multiplicity coordinates preserves that cone equivalence.

For the spin-half component, the raising matrix \(U\) has disjoint rows,
each with three coefficients \(s_i\in\{-1,1\}\), and \(UU^T=3I\).
The code constructs its full kernel without a learned dense rotation.
For each nonzero row, two integer kernel vectors are

\[
 (3s_1,-3s_0,0),\qquad (2s_0,2s_1,-4s_2).
\]

Zero columns contribute \(6e_i\). Dividing these vectors by six gives the
sparse rational basis \(K\). Exact integer checks establish
\(UK=0\), positive diagonal \(K^TK\), and
\(\operatorname{rank}K=\dim\ker U\). The spin-three-half highest weight
uses its full coordinate basis. All multiplicity cross-couplings remain.

The unreduced-magnetic control uses the same original word dictionary,
spin-averaged coefficient equations and zero initial Gram matrices. It has
the same attainable cone but a different parameterization and computational
cost. The separate first-64-columns control deliberately restricts that cone;
it is not claimed to select the best restricted basis.

Symmetry reduction of semidefinite and sum-of-squares programs is established
mathematics; see [Gatermann and Parrilo](https://arxiv.org/abs/math/0211450).
The present experiment tests a particular input-derived molecular realization,
its cold construction, cost and attained accuracy. It does not claim to invent
symmetry reduction.

## Exact upper and construction dependencies

The upper witness is an integer charge-conserving matrix product state.
Exact charge-flow checks establish its particle and spin-projection counts.
Integer tensor contraction computes

\[
 U=\frac{\langle\psi,H\psi\rangle}{\langle\psi,\psi\rangle}.
\]

No normalization approximation enters this rational quotient. The upper
witness does not need to be spin pure for the full fixed-N ground-energy
bound. Fresh runs discover this state using numerical DMRG, then construct
the nonsinglet proof, coefficient maps, spin projection, moment functional,
and singlet optimization afresh. No older MPS, Gram map, solution checkpoint,
factor certificate or FCI state is a constructor input. Stage 2 may restart
from stage 1 of that same fresh run; its cost is included.

The resulting guarantee is \(E_0\in[L,U]\). Numerical objective values,
residual norms and moment eigenvalues guide discovery only. A failed finite
search, a negative numerical dual eigenvalue, or a resource cap does not
establish an exact ceiling for the attainable certificate family.

## Sign-equivalent H8 fixtures

Canonical molecular orbitals can differ by signs between fresh integral
calculations. For a chosen set \(B\) of spatial orbitals, the unitary

\[
 V=\prod_{p\in B}(-1)^{n_{p\alpha}+n_{p\beta}}
\]

multiplies each fermionic word by the product of the corresponding orbital
signs. The bridge solves the sign constraints over \(\mathbb F_2\), then
checks the entire rational coefficient difference after transformation.
If that norm is \(\eta\), a new-fixture interval \([L,U]\) implies
\([L-\eta,U+\eta]\) for the old fixture. Here the H8 bridge finds orbital
indices 1 and 5 (zero based) and **\(\eta=0\)**: exact unitary equivalence.

## Fragment composition and its necessary sector condition

For mutually disjoint even number-conserving fragment Hamiltonians,
\(H=\sum_i H_i\), restricted to fixed local particle numbers \(N_i\),

\[
 \mathcal H=\bigotimes_i\mathcal H_{N_i},\qquad
 \sum_iL_i\le E_0(H)\le\sum_iU_i.
\]

The product of the fragment upper states has exactly the summed Rayleigh
energy. Lower operator inequalities tensor with identity and add. The
composition checker verifies disjoint orbital coverage, exact coefficient
addition, local-charge constraints, local charge boundaries of the product
MPS, and its full exact Rayleigh quotient. Added interfragment hopping and
removal of the local-charge conditions must be refused.

This proves interval composition for a decomposable problem. It does not
claim that a monolithic approximate family is size consistent, nor does it
handle charge redistribution by silently replacing the fixed local charges
with only a total charge. See also the distinct RDM size-consistency question
in [Nakata](https://arxiv.org/abs/1108.5665).

## Scaling and physical scope

Avoiding full determinant enumeration does not establish inexpensive scaling.
There are degree-six coefficient equations and cubic operator dictionaries;
their dimensions, sparse-map work, dense Gram eigensolves, state contractions,
exact integer sizes and convergence all matter. No uniform accuracy or
polynomial-time many-body solution theorem follows from these examples.

For physical differences, separately certified states give

\[
 E_S-E_T\in[L_S-U_T,\ U_S-L_T].
\]

The solver interval and physical-model uncertainty are separate quantities.
The accompanying CH2 geometry/basis study is a conventional numerical
CASSCF/SC-NEVPT2 control. It is retrospective, has no certified physical error
bar, and does not establish a better prospective experimental decision.
