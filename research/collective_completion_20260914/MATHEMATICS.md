# Exact scope of the collective construction

All accepted endpoints concern the supplied rational electronic Hamiltonian
on its entire fixed-particle-number sector. An MPS guides proposals and supplies
a separately verified upper bound. It never constrains the unknown ground state.

## Paired positivity and two-particle preparation

For odd cubic polynomials c_i, the sextic part of
c_i^dagger c_j + c_j c_i^dagger cancels by CAR. If S is real symmetric PSD,

    sum_ij S_ij (c_i^dagger c_j + c_j c_i^dagger)

is a sum of two positive quadratic forms and has degree at most four. Linear
words can be included. The sparse coefficient map can therefore be constructed
directly on quartic rows; its guiding moment matrix uses only the one- and
two-particle correlations. No three-particle moment table is necessary.

This is an established type of T1/T2 positivity, not a new general identity.
See [Mazziotti, 2012](https://arxiv.org/abs/1207.0541) and the
[v2RDM tutorial](https://arxiv.org/html/2310.10746v1) for the surrounding
representability framework. Our tested choice is a restricted, adaptive
congruence of these maps, with exact exported certificates.

The search uses S=V Q V^T with fixed rational V and Q PSD. Every retained
cross term is optimized subject to the Hamiltonian coefficient equations.
The numerical dual identifies negative directions of the omitted full moment
matrices. Appending their components outside V preserves the existing span.
Subspace exchange additionally truncates the numerical Q; this is only a
proposal and is not assumed to preserve a bound. Exact replay evaluates it.
The relevant precedent for subspace recovery is
[Ding et al., approximate complementarity](https://arxiv.org/abs/1902.03373).
The hypotheses of that paper are not asserted for these molecular instances.

## Independent linear responses

Consider the two polynomial lists (c,l_-,0) and (c^dagger,0,l_+), and apply
one common PSD matrix to both lists. Their cubic-cubic coefficients agree,
so their sextic terms cancel. The cubic-linear cross terms on the two sides
are independent. The resulting polynomial is still a positive sum of squares
of degree at most four. This is related to strengthened T2 constructions.

In discovery, its guiding moment matrix is partitioned as

    C = [[A,B],[B^T,D]].

A numerical pseudoinverse suggests the linear response -D^+ B^T and a
conditional cubic matrix A-B D^+ B^T. Low eigenvectors of that matrix select
cubic directions. The full linear response space is retained. This numerical
Schur calculation accepts nothing: exported integer polynomial factors are
checked against the exact original operator identity. Tests independently
check sextic cancellation, free quartic cross terms, and the full CAR map.

## Balanced-spin lifting

Let Hs be the exact rational SU(2) average of H, computed by the existing
adjoint-Casimir polynomial projector. Exact raising/lowering commutators
verify its invariance. Let delta=sum_w |(H-Hs)_w|. Every CAR monomial is a
contraction, so ||H-Hs|| <= delta on every sector.

For even N every allowed total spin is an integer. Each spin multiplet
contains M_S=0, and Hs is constant across the members of a multiplet.
Consequently its lowest fixed-N energy equals its lowest M_S=0 energy.
This does not assert that M_S=0 is a pure singlet sector.

Write Z=N_alpha-N/2. The checker verifies an ordinary global-N SOS proof
for Hs-ZY, where Y is Hermitian and preserves both spin counts. On M_S=0,
ZY vanishes. A lower bound b-eta for the auxiliary Hamiltonian therefore
gives the full original-H bound b-eta-delta. Eta is the exact coefficient
norm of the actual rational reconstruction residual. Tiny input spin defects
are paid, never silently dropped. The auxiliary Hamiltonian must match the
derived expression exactly; the checker refuses a substituted input.

## A singlet proof plus a collective nonsinglet proof

Every integer-spin multiplet with S>=1 has an M_S=1 member. A lower bound
Lt for Hs on M_S=1 therefore covers *all* nonsinglet spin multiplets.
Independently, on S=0 one has S^2=0, S_-=0, S_+=0, and M_S=0.

Thus the singlet checker may remove a scalar multiple of S^2 and an ideal

    S_+ W + W^dagger S_-,

where W is a number-conserving quadratic polynomial lowering M_S by one.
Its compression to the singlet subspace vanishes: P0 S_+=0 and S_- P0=0.
It need not vanish as an operator on the entire Fock space. The checker
requires the singlet label whenever these equalities are used.

If Ls and Lt are independently verified for Hs in the two pieces, then

    E0(H on N) >= min(Ls,Lt) - delta.

Both certificates are mandatory. Neither a guessed ground-state spin nor
an MPS spin expectation substitutes for the nonsinglet proof. The tiny
noninvariant perturbation in the original H is charged once at the end.

## Averaging the compact proof before reconstruction

For the screened singlet proof, apply the exact SU(2) projector T to the
whole proposed identity. T(SOS) is positive because it is a Haar average of
unitary conjugates of positive operators. The number ideal remains zero on
fixed N; the other admitted ideals have zero compression to the singlet,
which is invariant under the rotations. T(Hs)=Hs. Thus the verified residual
is T(R), and its exact coefficient norm is a valid allowance for the singlet
lower endpoint. The other-spin certificate and original-H delta remain
mandatory. A balanced M_S=0 proof alone cannot use this shortcut for arbitrary
spin ideals, because that subspace is not invariant under every rotation.

The coefficient projector is built from exact twirls of the independent
Hermitian quartic monomials. For H8 it is an integer matrix divided by six.
Its idempotence is checked exactly, with a proved integer overflow bound.
Its trace is the exact rational rank. The code selects that many independent
rows by elimination modulo the checked prime 65521. Modular independence
implies rational independence; idempotence and trace give the matching upper
rank bound. Consequently the 581 selected equations span all 1,525 projected
equations. This reduces equations without discarding any projected identity.

The final independent checker reconstructs the full polynomial and applies
the exact projector again. It does not trust the selected numerical rows or
their residual estimates. This is a compact proof *averaging recipe*; it does
not expand every rotated factor as a separate numerical optimization variable.

## Sharing spin orbits and orienting the numerical SDP

Each dictionary has a fixed particle charge and spin projection. Exchanging
alpha and beta differs from a pi spin rotation by a constant phase on such
a dictionary. That phase cancels in a square. Consequently

    T(p^dagger p) = T(flip(p)^dagger flip(p)).

For two dictionaries exchanged by this involution, their projected spans can
be placed in one representative dictionary. A numerical SVD proposes their
union; its truncation and rounding are not assumed exact. For a dictionary
fixed by the involution, separate positive and negative parity spans have
zero averaged cross term. The self-orbit splitting is checked algebraically.
All exported factors still undergo exact full-polynomial reconstruction.
Spin sharing is a reduction of redundant optimization variables, not a new
positivity theorem or an assertion that the trial state has exact spin.

The global-triple variant includes pure annihilation triples across the
whole orbital set, instead of confining that dictionary to the two fragments.
Their adjoints remain paired, so sextic cancellation and quartic preparation
are unchanged. No full many-body matrix is constructed.

The coefficient SDP can also be solved in the opposite orientation. Write
the SOS equations as F x + sum_k A_k(Q_k) = h, with Q_k positive and x_0=b.
Its moment dual is

    minimize h^T y, subject to F^T y=e_0 and A_k^*(y) positive.

Dual multipliers of that problem propose x and the Q_k. The off-diagonal
entries in the adjoint map carry a factor of one half when using one
coordinate per symmetric pair. A two-dimensional independently diagonalized
example checks this convention, coefficient reconstruction, and row scaling.
The moment problem has one variable per independent coefficient equation;
this changes numerical cost without changing the intended cone. Neither
primal-dual agreement nor numerical feasibility accepts an endpoint. The
ordinary integer-factor checker still reconstructs and pays the residual.

The independent-response search also supports an exact inclusion of a
previous paired span: copy its cubic coefficients into the common slots,
copy any tied linear coefficient into both response slots, and append the
independent linear coordinate vectors. This changes no old polynomial square.
It avoids confounding an additional response with replacement of the original
directions. The inclusion does not guarantee that a time-limited numerical
optimizer will recover an equally good certificate.

## Optional collective bound on reconstruction error

The actual residual can be bounded by the already implemented exterior-power
method instead of only its coefficient norm. After exact factor expansion
(and exact spin averaging where admitted), write its k-body part as a matrix
on k-orbital wedge indices. A rigorously checked lower eigenvalue ell_k gives
the fixed-N lower contribution binomial(N,k) ell_k. The accepted value is the
maximum of this bound and the coefficient-norm bound. Integer Gram witnesses
and rational Gershgorin remainders verify every small matrix inequality.

This witness is an additional proof dependency, with additional construction
and replay cost. Its k-orbital indices are not an enumeration of the N-electron
determinants. It only bounds the actual residual of an already reconstructed
SOS identity; it cannot legitimize an incorrect spin-sector or Hamiltonian
binding. Negative-energy examples and mutated wedge indices are tested.

## Coherent squares without a growing cubic Gram block

The atom variant keeps the full quadratic PSD blocks but writes each cubic
paired contribution as a nonnegative combination of fixed coherent squares:

    sum_j t_j (p_j^dagger p_j + p_j p_j^dagger),  t_j >= 0.

It can start from the two-particle MPS guide alone, or from the factors of an
already checked compact certificate. The latter has that certificate's entire
discovery as a dependency. Recorded power-of-two rescalings condition the
fixed directions. A missing block in a seed certificate contributes zero;
the implementation does not mistake that omission for a missing proof input.

At each iteration the coefficient dual is applied to the full paired moment
maps. A negative eigenvector proposes a new coherent polynomial p_j. This
tests an omitted positive-square inequality. It does not prove an exact
obstruction for the current family, and it does not guarantee a useful bound
improvement at each iteration. Existing atoms remain available. No full
determinant vector, Hamiltonian matrix, or previous full cubic Gram certificate
is used in this column-generation construction.

This trades a general PSD matrix on a fixed span for individually weighted
directions that can be added. The two finite cones differ. It avoids a growing
cubic PSD block, but adds scalar columns and may need many iterations. Factors
with identical word dictionaries are grouped only when serializing the final
proof; exact replay expands their actual squares and residual as before.

## Audit of the supplied spectral-filter theorem

The attached text's filter idea is a valid sufficient criterion. It has not
been implemented as an accepting molecular filter in this pass. The toy
archives referred to in that text were not available among the located files.
Polynomial eigenvalue-count filtering is established; e.g.
[Di Napoli, Polizzi and Saad](https://arxiv.org/abs/1308.4275) study stochastic
count estimates. Such estimates alone do not provide the deterministic bound
required here.

For completeness, write x_b=cosh(alpha)>1, x=cosh(beta)>=x_b, d_j=T_j(x_b).
On a hypothetical forbidden eigenvector, the stated symmetrized approximate
recurrence gives scalar errors r_j of magnitude <=epsilon_j and

    d_k g_k = T_k(x) + sum_{j=1}^k d_j U_{k-j}(x) r_j.

For j>=1,

    0 <= T_j(x_b) U_{k-j}(x)/T_k(x) <= coth(beta) <= coth(alpha).

One proof replaces T_j(x_b) by T_j(x), uses the hyperbolic product identities,
and bounds sinh((k-2j+1) beta) by sinh((k-1) beta). Since T_k(x)/d_k>=1,
if K sum epsilon_j<1, with K=coth(alpha), then

    g_k >= 1-K sum epsilon_j.

But |g_k| <= ||G_k P_N||_HS. Hence h+K sum epsilon_j<1 excludes the forbidden
eigenvalue. The norm is the **unnormalized** Hilbert-Schmidt norm. The exact
filter version counts each forbidden eigenvalue with weight at least one.
For any finite H and strict b<E0, choosing b<a<E0 and c>=Emax makes the exact
family eventually pass. This proves expressiveness, not cheap evaluation.

The code checks 1,376 scalar kernel cases with exact rational arithmetic.
On the molecular fixtures it reports only conditional degree/precision scales
using a crude certified spectral enclosure. These are not degree lower bounds,
not molecular trace certificates, and not a reason to rule out better filters.

## What the amplitude-ratio probe establishes

The molecular diagnostic queried 49 selected determinant labels and two
Hamiltonian columns per fixture. It did not enumerate either full sector.
Each fixture contains an exact positive-product triangle of off-diagonal
matrix elements. That cycle product is invariant under diagonal sign changes;
three negative edges would have negative product. Thus such a sign gauge
cannot make every transition negative in this occupation basis.

This does not obstruct general signed-Laplacian routing, another orbital
basis, or a phase-aware compact ansatz. The local spin-chain certificate in
the supplied text cannot simply be applied without constructing molecular
ratio rules and proving their collective capacity bounds. No new toy-chain
result is counted as molecular progress.
