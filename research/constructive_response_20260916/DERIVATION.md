# Verified mathematical results and their limits

This document is the consolidated, reviewed derivation for this pass. A proved
finite-dimensional inequality is not a proof of a tractable many-body class.
No priority or field-level breakthrough claim is made.

## 1. A sharp relative-residual energy bound

Let a finite Hermitian Hamiltonian, shifted by a real reference e, have blocks

\[
H-eI=\begin{pmatrix}A&B^*\\B&D\end{pmatrix},\qquad D\succeq\delta I,
\quad\delta>0.
\]

For any trial response X from the retained to the eliminated space, define

\[
R=B-DX,\quad M=I+X^*X,\quad
K=A-B^*X-X^*B+X^*DX.
\]

Suppose two operator inequalities have been certified:

\[
K\succeq0,\qquad R^*R\preceq\rho^2 M,
\qquad 0\le2\rho<\delta.
\]

Then

\[
\boxed{H\succeq(e-\eta)I,\qquad
\eta=\frac{\rho^2}{\delta-2\rho}.}
\]

Proof: the invertible triangular map T=[[I,0],[-X,I]] gives

\[
T^*(H-eI+\eta I)T=
\begin{pmatrix}
K+\eta M&R^*-\eta X^*\\R-\eta X&D+\eta I
\end{pmatrix}.
\]

For every retained vector v,
||Rv||<=rho sqrt(v*Mv) and ||Xv||<=sqrt(v*Mv). Hence

\[
(R-\eta X)^*(D+\eta I)^{-1}(R-\eta X)
\preceq\frac{(\rho+\eta)^2}{\delta+\eta}M=\eta M.
\]

The last equality is rho^2=eta(delta-2rho). The Schur complement is PSD,
so the full congruence and therefore H-e+eta are PSD. Nothing commutes by
assumption. The only inverse is in the proof; the checker can verify the gap
and the two displayed quadratic inequalities without forming that inverse.

The metric M is exactly the norm of the lifted vector [v;-Xv]. Thus the bound
does not directly penalize a large ||X||. It still requires a *global operator*
bound on the residual in that metric. Sampling a few v is insufficient.

### The constant and threshold are optimal for this information

For scalar retained and eliminated spaces choose, for arbitrarily large k,

\[
D=\delta,\ X=k,\ R=-\rho k,\ B=(\delta-\rho)k,
\ A=(\delta-2\rho)k^2.
\]

Then K=0 and R^2<=rho^2(1+k^2). For any proposed nonnegative allowance h,

\[
\det(H-eI+hI)
=k^2[(\delta-2\rho)h-\rho^2]+h(\delta+h).
\]

If delta>2rho and h<rho^2/(delta-2rho), this determinant is negative for
sufficiently large k. Thus no smaller uniform allowance follows from these
same hypotheses. If delta<=2rho (with rho>0), every finite h can be defeated.
`sharpness_counterexample` chooses an integer k by exact arithmetic and checks
the original premises and negative determinant independently.

This is an information limit of the abstract certificate. A specified physical
Hamiltonian supplies other information; its actual ground energy can be better
bounded by using that information. Failure of this rule is not physical
intractability or a limit of all response methods.

### Exact upper bound and the numerical control

A nonzero retained trial v gives a variational upper through [v;-Xv]:

\[
U=e+\frac{v^*Kv}{v^*Mv}.
\]

The four-site square test has U/t=8, four hopping edges of strength t=1,
N_up=N_down=2, and a retained space with one electron per site. Each local bond
has hopping >=-2I on **all** its occupation sectors. Because the eliminated
space has at least one doublon,

\[
QHQ\succeq(8-4\times2)Q=0,
\quad D=Q(H-e)Q\succeq(-e)Q\quad(e<0).
\]

The gap proof uses a 16-dimensional local bond check, independently constructed
from the CAR. It does not diagonalize QHQ. The test nevertheless enumerates all
36 balanced global configurations to construct the response and its 6-by-6
retained Gram matrices. It is a correctness control, not the missing compact
constructor. Cold discovery uses 16 Jacobi updates of a 30-by-6 response and
two six-dimensional generalized Ritz solves, without a full ground-state
eigensolver. Rounded proposals are checked with rational arithmetic.

The accepted interval is approximately

\[
[-1.320235241487274,\;-1.320234922543856]t,
\quad\text{width }3.1894341811\times10^{-7}t.
\]

The exact rational endpoints are in the replay receipt. Here t is the Hubbard
hopping energy; these numbers are not molecular hartrees or a chemical-accuracy
claim. The relative residual allowance is 2.4148727373e-7 t versus an independently
checked absolute-residual allowance of 2.7678424674e-7 t. This small example
supports the implementation, not a competitive-advantage claim.

The proof uses established Schur, congruence, and residual methods. Rigorous
block resolvent enclosures already exist, for example
[Zimmerling, Druskin and Simoncini](https://doi.org/10.1007/s10915-025-02799-z).
The exact normalization and sharpness above were derived in this pass;
literature priority is unresolved.

## 2. A global bare projector can become poorly conditioned

For a normalized ground vector psi, set u=Ppsi, v=Qpsi, p=||u||^2, C=QHQ.
For 0<p<1, the projected eigenvalue equations give

\[
v^*(C-E_0)v=u^*(PHP-E_0)u,
\quad
\lambda_{\min}(C)-E_0\le\frac{u^*(PHP-E_0)u}{1-p}.
\]

If the ground state is unique with full-system gap g>0, every normalized Q
vector has ground overlap squared at most 1-p. Consequently

\[
\lambda_{\min}(C)-E_0\ge gp>0.
\]

The inverse at E0 therefore exists, and the projected equation yields

\[
\|(C-E_0)^{-1}QHP\|\ge\sqrt{(1-p)/p}.
\]

For L independent Hubbard dimers at U/t=8,

\[
e_d=4-2\sqrt5,\quad w=(1+2/\sqrt5)/2,\quad E_0=Le_d,\ p=w^L.
\]

The unique product ground and its gap g=-e_d remain valid with only global
charge and spin projection fixed. To see why charge transfer does not create
a lower state, subtract the chemical-potential term 4(N_dimer-2) from each
local dimer. The sector minima for local particle numbers 0,1,2,3,4 are
8,3,e_d,3,8. The next two-electron energy is zero, realized by a triplet.
The chemical-potential shifts cancel at total N=2L. Thus the full balanced
sector has unique ground and gap -e_d; it is not a product-sector restriction.

For the global bare no-doublon projector PHP=0, so

\[
(-e_d)w^L\le\lambda_{\min}(C)-E_0
\le\frac{(-Le_d)w^L}{1-w^L}.
\]

Both bounds decay exponentially up to the displayed linear factor. At L=256,
the upper bound is about 1.1297327354e-4 t and the response norm is at least
1034.3455. The full system gap remains about 0.472136 t. This disproves a
uniform-conditioning inference for that global bare projector, even on a
physically easy product family. It does not establish exponential complexity
at a fixed energy tolerance: moving below E0 keeps a resolvent separation, and
energy-weighted bounds can differ substantially from absolute response bounds.

The exactly rational local dressing used in the probes is specified by
r=116434/10^6, c=(1-r^2)/(1+r^2), s=2r/(1+r^2). For one dimer, the overlap
with the exact ground is

\[
w_d=\tfrac12+
\frac{4(c^2-s^2)+4cs}{\sqrt{80}},
\quad a_d=-4cs+8s^2.
\]

Its product gap *upper* estimate must use
L(a_d-e_d)w_d^L/(1-w_d^L), not the bare expression. The local rotation is
almost, but not exactly, the ground-state rotation. All radical and probability
intervals in the replay are rational and outward rounded.

## 3. Composition and local coercivity do not close the hard part

For an ordered full-operator enclosure M_-<=M<=M_+, Schur minimization
preserves that order when the corresponding eliminated blocks are positive.
Two successive exact eliminations equal elimination of their union. The
five-dimensional exact regression checks this with perturbations in diagonal
and cross blocks. Error sizes need not merely add: the positive matrices
[[40,6],[6,1]] and [[40,6],[6,2]] have Schur complements 4 and 22. A unit change
in the eliminated block caused an eighteen-unit reduced change.

For one Hubbard bond the optimal scalar in T_ij+a(d_i+d_j)+bI>=0 is

\[
b(a)=\max\{t,(\sqrt{a^2+16t^2}-a)/2\},\quad a,t\ge0.
\]

The one- and three-electron sectors must be included. At degree z and edge
count m, summing this inequality gives
H>=(U-az)D-mbI; replacing D by I on Q requires U-az>=0. This is a valid
certificate but its extensive loss is generally too large. The successful
four-site gap above does not imply a gap rule for general 2D lattices.

The still unproved result would construct X, a retained positivity proof for K,
and the residual-metric inequality directly in a representation whose cost and
verification remain controlled under interacting reductions. None of the
identities or finite checks in this pass establishes that result.
