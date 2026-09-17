# General many-body objective and exact sector-action quotient

## Objective and status

The research objective remains a general solution strategy for interacting
quantum many-body systems. A sequence of increasing hydrogen chains is not that
objective. The original question remains: can the difficult correlations be
discovered as a manageable collection of interference patterns, while everything
else is eliminated or certified collectively?

There are distinct mathematical obligations: energy and spectral guarantees,
state or observable reconstruction, real-time evolution, and statistical
mechanics. A ground-energy interval alone does not discharge the others.
Infinite-dimensional and continuum problems also require controlled truncation
or separate analysis; solving a finite orbital Hamiltonian does not settle them.

The substantive unresolved requirement is a **constructible representation**
with controlled errors and resource growth. In-principle completeness with
exponential resources already exists. A universal efficient method for arbitrary
local Hamiltonians would have major complexity consequences; known QMA results
are conditional barriers, not an unconditional impossibility theorem. See the
[primary sources](SOURCES.md). Neither a universal solution nor a general scaling
theorem is claimed here.

This pass implements one reusable exact algebraic component. Its theorem is a
faithful-trace/linear-algebra argument, not a claimed new solution of the
many-body problem. It removes *exactly redundant* operator directions. It does
not identify all the *energetically important* nonredundant ones.

## 1. The complete action kernel of any finite dictionary

Let P be the orthogonal projector onto a nonempty, finite-dimensional physical
sector, with d = Tr P. For any operator A, define

\[
\tau_P(A)=d^{-1}\operatorname{Tr}(PAP).
\]

For a finite operator dictionary O_1,...,O_m, form

\[
G_{ij}=\tau_P(O_i^\dagger O_j).
\]

For any coefficient vector z, set Z = sum_i z_i O_i. Then

\[
z^\dagger Gz=\tau_P(Z^\dagger Z)
=d^{-1}\|ZP\|_F^2.
\]

Consequently G is positive semidefinite and

\[
\boxed{\ker G=\{z:(\sum_i z_iO_i)P=0\}.}
\]

This is an equality, not a sufficient subset of null directions. In particular,
it detects operators vanishing on every physical state in the sector, including
states far from the trial MPS. No energy or particular ground state appears.

Select pivot operators W_a = O_{p_a} from an exact positive LDL decomposition of
G. They span the quotient by ker G. Exact coordinates C obey

\[
OP=C^TWP,\qquad G=C^\dagger G_{pp}C,\qquad G_{pp}>0,
\]

where O and W are column lists; for complex coefficients the first expression
uses the ordinary transpose and the Gram identity the conjugate transpose. The
implementation uses real rational coefficients only. Every difference
O_j - sum_a C_{aj}W_a has zero action on the sector. Completeness is certified by
the exact positive pivots and one independent kernel vector per free column.

## 2. Compute the singlet trace without listing the sector

For s spatial orbitals and N=2n particles, write T for exact averaging over spin
rotations. The normalized singlet trace is

\[
\tau_0(A)=\frac{\operatorname{Tr}_{N,M_S=0}T(A)
-\operatorname{Tr}_{N,M_S=1}T(A)}{d_0},\qquad
d_0={s\choose n}^2-{s\choose n+1}{s\choose n-1}.
\]

Spin averaging is essential for general products. In each integer-spin
multiplet, an invariant operator has the same trace in each magnetic component;
subtracting the two components cancels all S>=1 multiplets. Only S=0 remains.

For a normal-ordered monomial, a magnetic trace vanishes unless its creation
and annihilation index sets agree. For matching sets containing a alpha and b
beta orbitals, occupation counting gives

\[
\operatorname{Tr}_{N_\alpha,N_\beta} E
=\sigma_E {s-a\choose N_\alpha-a}{s-b\choose N_\beta-b},
\]

where sigma_E is the sign required by the actual annihilator-order convention.
The existing CAR implementation stores both creation and annihilation indices
in increasing order, giving sigma_E = (-1)^{k(k-1)/2} for k=a+b.

The implementation reuses the existing exact degree-six spin twirl and trace;
therefore its input dictionaries currently have degree at most three. The
general action-kernel theorem does not impose this degree limit. Extending the
implementation to higher degree requires additional trace machinery and cost
accounting, not merely changing a validator.

## 3. Additional exact sparsity, with no locality approximation

Each supplied operator must have a single spatial charge vector q: conjugation
by independently chosen spatial-orbital phases sends O to exp(i q.theta) O.
These rotations commute with particle number and total spin, and preserve the
uniform singlet trace. Hence

\[
q_i\ne q_j\quad\Longrightarrow\quad\tau_0(O_i^\dagger O_j)=0.
\]

This partitions the trace Gram matrix exactly. It does **not** delete couplings
from a Hamiltonian or assert that the molecular Hamiltonian has these phase
symmetries. The symmetry belongs to the reference trace.

Only within-charge products are generated. If the groups have sizes m_q, the
number of trace products is sum_q m_q(m_q+1)/2 rather than m(m+1)/2. The algorithm
still pays for those products, rational elimination, coefficient bit lengths,
and the original dictionary. These counts are not whole-solver runtime claims.

## 4. What can and cannot be removed from an SOS problem

For the real coordinates C used here and any real PSD Gram coefficient Q,

\[
P O^\dagger Q OP=P W^\dagger(C Q C^\dagger)WP,
\qquad CQC^\dagger\succeq0
\]

For complex coordinates the reduced matrix is instead conjugate(C) Q C^T,
which is also PSD. Thus the quotient preserves the physical
positive-square cone, **modulo all the displayed sector-annihilating cross
identities**. Cross terms O_i^dagger Z + Z^dagger O_i vanish between P's whenever
ZP=0.

That is not automatically equivalence to an old truncated coefficient SDP. An
old SDP may omit some of those identities. To claim equivalence one must prove
that the required cross identities lie in its admitted ideal span, excluding
the free energy shift. Alternatively one may add validated sector identities,
declare a strengthened relaxation, and extend the accepting checker explicitly.
This pass does neither; existing molecular acceptance remains unchanged.

A simple counterexample explains the gate. Take P=|0><0|, O=(I,Z), with Z the
Pauli matrix. The physical trace Gram is [[1,1],[1,1]], so (1,-1) is null. But
the full-space normalized functional Tr(A)/2 has Gram I_2. It obeys positivity
and Z^2=I while failing the sector relation ZP=P. A physical null relation is
not forced by every relaxation that omitted the relation defining the sector.

## 5. Consequence for exact dual repair

Once the appropriate relations are admitted and imposed, the trace is strictly
positive on the reduced dictionary. If a candidate functional y already obeys
all other required affine and residual inequalities exactly, a rational alpha
with

\[
M_y+\alpha G_{pp}\succeq0
\]

can be checked by exact PSD arithmetic. The convex combination

\[
y'=(y+\alpha\tau_0)/(1+\alpha)
\]

then preserves the common feasible affine/convex constraints, provided tau_0
also satisfies those constraints. One must use a common alpha valid for all
blocks. There is no claim that this mixing fixes an uncorrected affine equation
or an inequality on whose boundary tau_0 lies. Its energy changes by

\[
y'(H)-y(H)=\frac{\alpha}{1+\alpha}[\tau_0(H)-y(H)].
\]

That cost can destroy a useful obstruction margin. Exact null discovery solves
one prerequisite, not the complete dual repair.

## 6. The larger missing construction remains explicit

The earlier response work already establishes, for a retained/eliminated split
of H-bI with blocks A,B,D and D>=delta I, that an approximate response X gives

\[
K=A-BX-X^\dagger B^\dagger+X^\dagger DX,\quad
R=B^\dagger-DX,\quad
S=K-R^\dagger D^{-1}R.
\]

A controlled residual gives S>=K-eta I. The unresolved constructive obligation
is to discover and certify K-eta I>=0 economically, then keep the induced
operators manageable under further elimination. Borrowing a successful global
SOS proof to certify K transfers a proof but does not solve that discovery
problem. See the [existing response derivation](../global_response_20260913/DERIVATION.md).

The action quotient above can remove algebraic redundancy before that search.
It does not prove that the remaining correlations admit a small dictionary,
that such a dictionary can be found cheaply, or that the representation closes
under energy elimination, observables, and dynamics. Those are the central
general-mechanism questions. They will be evaluated as mathematical and
computational obligations, not replaced by a requirement to reach the next H
number.
