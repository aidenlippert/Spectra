# Exact joint-pattern certificates and their obstruction

All statements concern the frozen rational electronic H6 Hamiltonian, twelve
spin orbitals, and the full fixed-N=6 sector. The H4 control uses N=4 in eight
spin orbitals. No occupation-sector matrix enters discovery or acceptance.

## The retained Hamiltonian is fixed

The preceding molecular pass gives ten real Hermitian density operators

\[
Q_k=\sum_{pq,\sigma}L^k_{pq}a^\dagger_{p\sigma}a_{q\sigma},\qquad
H_{10}=\widetilde h+\tfrac12\sum_{k=1}^{10}\lambda_k Q_k^2.
\]

Its rational tail certificate proves, on N=6,

\[
H_{10}+\ell I\preceq H\preceq H_{10}+uI,\qquad
\ell=-3.632\times10^{-9},\quad u=1.08327632\times10^{-4}
\]

in Hartree. Both the original fixture and tail are hash-bound into every new
certificate. The accepting path replays the coefficient PSD checks and number
identity from [the preceding proof](/Users/aidenlippert/Documents/Spectra/research/molecular_collective_20260913/PROOF.md).
Even the two- and four-pattern constraint trials use this same H10.

## Joint fermionic constraints

For a chosen subset of patterns, use the odd generators a_i and Q_k a_i.
The latter are normal-ordered cubic polynomials. Within each exact spin-charge
and orbital-parity symmetry class, let B be any real linear combination of
all these generators. Then

\[
\{B^\dagger,B\}=B^\dagger B+BB^\dagger\succeq0.
\]

For generators P_i, the corresponding matrix-valued polynomial is
S_ij=P_i^dagger P_j+P_j P_i^dagger. A real PSD Gram matrix produces a sum of
these positive anticommutators. The map for an off-diagonal Gram entry includes
S_ij+S_ji; forgetting this factor of two would change the optimization.

The degree-six terms cancel exactly. At leading degree, interchanging two
cubic odd monomials gives sign (-1)^9=-1. CAR contractions lower the degree.
Consequently their anticommutator has degree at most four. Terms involving
linear generators have no higher degree. The implementation computes both
ordered products and checks the cancellation before forming the numerical map.
It never discards a surviving degree-six coefficient.

The quadratic baseline consists of all linear, particle-particle, hole-hole,
and particle-hole generators, partitioned by exact molecular symmetries.
The joint H6 addition has four blocks of dimension 33: three linear generators
and thirty pattern-times-annihilator generators per block. The separate-pattern
ablation embeds each pattern and the linears in its own block, hence permits
no cross-pattern Gram entries. This is a subcone of the joint construction.

## Exact accepting lower bound

A compact certificate expands into the identity

\[
H_{10}-bI=\sum_j T_j^\dagger T_j
 +\sum_r\{B_r^\dagger,B_r\}+(\widehat N-6)X+R,
\]

where X is Hermitian and at most one-body, and T_j are at most quadratic.
All coefficients and factors are rational. On the chosen sector the number
ideal vanishes. Every normal-ordered CAR monomial has operator norm at most
one, so if R=sum_w r_w w,

\[
E_0(H_{10})\ge b-\sum_w|r_w|,\qquad
E_0(H)\ge b-\sum_w|r_w|+\ell.
\]

Numerical Gram matrices are PSD-clipped and factored only as proposals.
For each anticommutator, one rounded factor row defines B; its partner B^dagger
is obtained by exact adjunction. Both are sent to the unchanged rational SOS
verifier. Thus rounding cannot destroy the degree-six cancellation. All
remaining coefficient error, including clipping and solve error, is charged
in R. An `optimal_inaccurate` solver status is not an acceptance condition.

## Exact dual ceiling

Let y be a real Hermitian linear functional on balanced CAR monomials through
degree four, with y(I)=1 and |y(w)|<=1. Require y((N-6)X)=0 for every one-body
Hermitian X, and require positivity of all quadratic and joint-pattern moment
Gram matrices. Then for any certificate in this cone,

\[
y(H_{10})=b+y(\text{positive terms})+y(R)
\ge b-\sum_w|r_w|.
\]

Thus y(H10) is a ceiling on the entire lower-bound family, independently of
whether a numerical solver has converged. No strong-duality assumption is used.
The witness is zero on all nontrivial symmetry classes. Rebuilding the complete
generators therefore gives block-diagonal moment matrices with exactly the
checked blocks; off-block moments are zero. All one-body number-ideal identities
are checked, including those outside the numerical symmetry support.

The saved numerical proposal is rounded, made Hermitian, and projected onto the
number identities by exact affine elimination. A mixture of 10^-7 with the
uniform fixed-N trace functional is accepted. Smaller tested mixtures fail
an exact PSD pivot test. For a canonical diagonal k-body word, that trace is
(-1)^(k(k-1)/2) binomial(N,k)/binomial(m,k); off-diagonal moments vanish. This
seed requires no enumeration of determinants.

The accepted witness has 919 moments. It passes 26 quadratic PSD checks,
four 33-dimensional anticommutator PSD checks, normalization, the coefficient
box, and 79 exact number-ideal identities. Integer coefficient congruences
produce its Gram matrices; rational symmetric elimination, including the
zero-pivot row condition, decides PSD.

Its molecular lower-bound ceiling is

\[
C=y(H_{10})+\ell=-6.346068751063810\ldots\ \mathrm{Ha}.
\]

The independently replayed reference upper is
U=-6.333058626233001... Ha. Therefore every interval using this cone and that
upper has width at least U-C=13.0101248308097... mHa. The best new lower lies
only 0.0003327866847... mHa below C. The ceiling concerns the full ten-pattern
cone, which includes the smaller joint two-pattern construction by embedding
its Gram matrices and setting other coefficients to zero.

There is also a physical error obstruction. The prior full cubic precision
certificate, including its separate spectral residual proof, is replayed
against the exact same H and N. It proves

\[
E_0(H)\ge L_*=-6.3331136206258359\ldots\ \mathrm{Ha}.
\]

For every lower bound L produced by the present family,

\[
E_0(H)-L\ge L_*-C
=12.9551304379743\ldots\ \mathrm{mHa}>1.6\ \mathrm{mHa}.
\]

This last argument uses the old physical proof only as a comparison witness.
Neither it nor its discovery factors enter the new pattern constraints or SDP.
Its construction cost and the prior upper-witness discovery remain separate.

## A concrete missing spin constraint

Split each existing spatial pattern into its two spin components:
Q_k=Q_k,alpha+Q_k,beta. No new spatial coefficient matrix is introduced.
Use a_i and Q_k,sigma a_i to construct four diagnostic blocks of dimension 63.
Their most negative numerical eigenvectors are proposals only.

The exported integer direction z gives a polynomial B and passes the exact
test

\[
\frac{y(\{B^\dagger,B\})}{\sum_i z_i^2}
=-0.0115165170755110\ldots<0.
\]

This number is a normalized moment diagnostic, not an energy improvement.
The direction has 63 nonzero coefficients and 156 operator monomials; its
anticommutator has 919 monomials, all of degree at most four. Its expectation
would be nonnegative in every physical state, so it separates the accepted
dual functional from a necessary fermionic positivity condition.

The spin-resolved cone was not optimized. Another dual functional might retain
the same low objective after this particular witness is removed. A new bound
calculation is necessary to establish any energy gain.

The obstruction excludes only the defined quadratic-plus-pattern
anticommutator cone with the stated number ideal and coefficient-L1 residual
bound. It does not exclude spin-resolved constraints, other correlated
generators, different residual certificates, or compact many-body methods.
