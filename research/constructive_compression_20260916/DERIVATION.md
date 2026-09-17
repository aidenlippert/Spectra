# Constructive tensor certificate: proved scope and unresolved step

This pass establishes a checkable spectral-witness construction. It does not
establish efficient compression of general many-body systems. The target
0.001t interval on the fully coupled eight-site Hubbard ladder was not reached
by the new lower construction. Previously accepted certificates are preserved.

## 1. Declared physical model and exact seed

The open 2×r ladder has 2r sites, 3r−2 edges, hopping t=1, U=8 and
N_up=N_down=r. Apply the bipartite particle–hole transformation on down modes,
with eta_(a,s)=(-1)^(a+s). In the transformed balanced sector,

\[
 H'=K_\uparrow+K_h-8\sum_i n_{i\uparrow}n_{ih}+8rI.
\]

The constant is 8N_up before restricting the sector. Replacing it by 8r on
all of Fock space would be wrong. The validation checks this transformation on
the fixed sector, including the CAR signs.

Regroup the fermions into up and hole configurations. With
\(d=\binom{2r}{r}\), vectors become d×d coefficient matrices and

\[
 \mathcal L(C)=KC+CK-8\sum_iN_iCN_i+8rC.
\]

K is real symmetric and N_i are real diagonal occupation projectors. The
map is self-adjoint for the Hilbert–Schmidt inner product. It has a normalized
positive-semidefinite ground eigenmatrix W. Here is the relevant variational
argument, without a claim of novelty:

* Decompose a general C as X+iY with X,Y Hermitian. Its norm and quadratic
  energy split into the sums for X and Y, so a Hermitian minimizer exists.
* Replacing Hermitian C by |C| preserves C² and the norm and kinetic energy.
* In C's eigenbasis,
  \(\operatorname{Tr}(CN_iCN_i)=\sum_{ab}\lambda_a\lambda_b|(N_i)_{ab}|^2
  \le\operatorname{Tr}(|C|N_i|C|N_i)\).
  The coefficient −8 therefore makes the energy no larger.

Consequently choose W≥0, ||W||_HS=1 and L(W)=E0 W. The unnormalized identity
seed obeys

\[
 \langle W,I_d\rangle=\operatorname{Tr}W\ge1,
 \qquad\|I_d\|_{HS}=\sqrt d.
\]

This is the spin-reflection-positive setting of
[Lieb's Hubbard theorem](https://doi.org/10.1103/PhysRevLett.62.1201), reviewed
with its assumptions by
[Boretsky, Cohn and Freericks](https://arxiv.org/abs/1712.02694).
The implementation represents I_d by a charge counter: each spatial site is
00 or 11 in up/hole occupation, and exactly r sites are occupied. Its bond
labels count occupations. No configuration list is constructed. Regrouping
introduces the constant phase (-1)^(r(r−1)/2), which is irrelevant to the
absolute overlap; for r=4 it is +1.

## 2. Inexact Chebyshev lower-energy theorem

Let H be any finite Hermitian operator with a normalized ground vector g.
Assume an explicitly justified seed overlap |<g,v0>|≥gamma>0. Choose real
b>ell and rational 0<z<1, and define

\[
 A=\frac{b+\ell-2H}{b-\ell},\qquad
 L=\frac{b+\ell}{2}-\frac{b-\ell}{4}(z+z^{-1}).
\]

b does not need to be a spectral upper bound for the logic of this theorem.
It affects filter stability and cost. The experiment uses the safe operator
bound b=8r+2(3r−2), equal to 52 for r=4.

The stored vectors are arbitrary exact rational tensors. Define their actual
residuals, rather than trusting a numerical compression routine:

\[
 r_1=v_1-Av_0,\qquad
 r_j=v_j-2Av_{j-1}+v_{j-2}\quad(j\ge2).
\]

Obtain rigorous upper bounds ||r_j||≤eta_j and ||v_k||≤nu_k. Then

\[
 \boxed{
 2z^k\nu_k+\frac{2}{1-z^2}\sum_{j=1}^kz^j\eta_j<\gamma
 \quad\Longrightarrow\quad E_0\ge L.
 }
\]

Proof: if E0<L, its A-eigenvalue is
x0>(z+z^(-1))/2>1. Projecting the recurrence onto g gives

\[
 T_k(x_0)\langle g,v_0\rangle
 =\langle g,v_k\rangle-\sum_{j=1}^k
 U_{k-j}(x_0)\langle g,r_j\rangle.
\]

Write x0=cosh(theta0) and z0=exp(−theta0)<z. The exact hyperbolic formulas give

\[
 T_k(x_0)^{-1}\le2z^k,\qquad
 \frac{U_{k-j}(x_0)}{T_k(x_0)}
 \le\frac{2z^j}{1-z^2}.
\]

Cauchy–Schwarz contradicts the boxed strict inequality. A seed approximation
with certified error eta0 adds eta0 to the left side. The implemented
Hubbard checker requires the exact seed and uses gamma=1, eta0=0.

An optional sharper finite-k theorem replaces the coefficients by

\[
 a_k=\frac{2z^k}{1+z^{2k}},\qquad
 w_{kj}=\frac{2z^j(1-z^{2(k-j+1)})}
 {(1-z^2)(1+z^{2k})}.
\]

These bounds are also valid at the unknown x0. To see the required ratio
monotonicity, put m=k−j+1≤k. The logarithmic derivative of
U_(k−j)(cosh theta)/T_k(cosh theta) is
m coth(m theta)−coth(theta)−k tanh(k theta). Since q coth(q theta)
increases with q and sinh(2k theta)≥k sinh(2theta), it is at most
2k/sinh(2k theta)−coth(theta)≤−tanh(theta)<0.
The implemented gate deliberately uses the simpler geometric envelope. A
failure of that envelope is not a failure theorem for the sharper gate.

Reducing b is not uniformly beneficial: it can make unwanted excited
components grow outside the filter's interval and increase ||v_k||.

## 3. Exact tensor replay

Every local tensor entry is an integer divided by a declared denominator.
The checker validates the charge flow and physical-model hash. It then:

1. Constructs H v_(j−1) from the local Hamiltonian MPO.
2. Contracts all pairwise inner products in the recurrence residual using
   Python integers and rational arithmetic.
3. Forms the exact squared residual, including cancellation between terms.
4. Uses integer square root and an outward rational ceiling for eta_j.
5. Evaluates the scalar gate using fractions only.

No global configuration list, dense many-body vector, dense Hamiltonian or
PSD diagonalization occurs on this path. The SVD outputs are untrusted input
proposals. Construction, norm reconstruction, rational bit length and tensor
transfer storage all still cost work and are recorded.

The early-refusal rule is narrowly defined. Once the geometric-envelope
error term alone reaches gamma, extending the same recurrence prefix cannot
make that specific gate pass at the same b,ell,z. This does not prove that a
different compression, seed, scaling, or proof rule must fail.

## 4. Preparation does not remove the error obligation

Given a certified ground-energy upper U<b, let B=(b−H)/(b−U). Its ground
eigenvalue is at least one. If
u_j=a_j B u_(j−1)+r_j, a_j>0, and a real ground projection of u_(j−1) is at
least gamma_(j−1), a valid next lower overlap is

\[
 \gamma_j=a_j\gamma_{j-1}-\|r_j\|.
\]

Normalizing a proposed compressed B u by n gives
(gamma_previous−eta)/n. This must remain positive. The numerical warmup
screens track this formula but their SVD error estimates are not certificates.

Orthogonality of a discarded SVD component does not justify quadratic error
accumulation against an unknown ground vector. Also, u=q+s with s≥0 implies
<W,q>≤<W,u>, so discarding positive mass does not preserve a lower overlap.
An appropriate certified operator ordering or additional spectral information
would be needed to improve these guarantees.

## 5. Mathematical extension and its limit

For arbitrary Hermitian H on a finite sector of dimension D, work in operator
space with L(X)=HX and Hilbert–Schmidt norm. L is self-adjoint, and for any
normalized ground g the rank-one projector W=|g><g| satisfies
L(W)=E0 W, ||W||_HS=1 and <W,I_sector>=1. Thus the same theorem applies with
an identity-operator seed, without Hubbard reflection positivity.

This mathematical extension is not implemented by the present real,
physical-dimension-two Hubbard checker. It generally requires an operator
MPO and physical dimension squared. The normalized seed overlap is 1/sqrt(D).
For an exact filter with bounded excited components, that overlap affects
degree logarithmically, but the needed MPO bonds and certified truncation
accuracy can grow very badly. No bound on those costs has been proved here.

The unresolved breakthrough is a constructor that keeps the weighted residual
small at useful precision while keeping representation and complete cost
manageable on growing interacting systems. The theorem is a way to check
such a constructor; it does not establish that one exists for general H.
