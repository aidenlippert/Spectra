# A perturbative order-versus-error guarantee

For the matched reference, full defect-Dicke enrichment admits an explicit accuracy bound. This theorem concerns **full enrichment at each order**, not the much smaller single-direction targeting heuristic. A separate [residual-targeting theorem](marginal_targeting_convergence.md) now supplies conditional selective convergence using an excited-state gap and explicit approximate-Ritz errors.

Let H=H0+delta, let P be an H0-invariant reference subspace, and write Q=I-P. Assume certified bounds

\[
QH_0Q\succeq cQ,\quad\|\delta\|\le\eta,\quad
b\le E_0(H)\le u<c-\eta,\quad\lambda_{\max}(H)\le B.
\]

Set rho=eta/(c-u)<1 and W=B-b. Define U0=P and let Ur be the H0-invariant closure of U(r-1)+delta U(r-1). If Er is the lowest Rayleigh value in Ur, then

\[
\boxed{0\le E_r-E_0(H)\le W\rho^{2r+2}.}
\]

This is a sufficient bound. It does not claim that the order it prescribes is necessary or that the full enriched space is small in practice.

## Proof

Take a normalized ground vector psi with energy E. The complementary resolvent R=(QH0Q-E)^-1 exists. Since E<=u,

\[
\|RQ\delta Q\|\le\rho,\qquad
Q\psi=-RQ\delta\psi,\qquad\|Q\psi\|\le\rho.
\]

Writing A=-RQ delta Q and s=-RQ delta P psi, repeated substitution gives

\[
Q\psi=\sum_{j=0}^{r-1}A^js+A^rQ\psi.
\]

Each retained term contains at most r perturbation factors and lies in Ur. Projection by Q only subtracts a vector in P. The resolvent preserves every H0-invariant extension of P inside Q; in finite dimension it acts within the corresponding Krylov closure. Therefore

\[
\operatorname{dist}(\psi,U_r)\le\|A^rQ\psi\|
\le\rho^{r+1}.
\]

Let e be the true distance and normalize the orthogonal projection of psi onto Ur to obtain phi. Then |<psi,phi>|²=1-e². The component of phi orthogonal to psi has squared norm e², so

\[
0\le\langle\phi,H\phi\rangle-E
\le(\lambda_{\max}(H)-E)e^2
\le W\rho^{2r+2}.
\]

Minimizing over Ur gives the claim. No commutation between delta and H0 is used. Ground-state degeneracy does not invalidate the argument.

## Connection to the implicit dictionary

Each degree-at-most-four perturbation monomial touches at most four pair labels. H0 evolution never enlarges the distinguished set of a defect-Dicke atom. Thus full order r is representable in the dictionary

\[
\sum_{q=0}^{\min(4r,m)}\binom mq4^q(m-q+1).
\]

The dictionary is polynomial in m for fixed r. The error theorem now supplies a sufficient r when rho and W are controlled. It does **not** make the resulting cost uniformly polynomial in system size and inverse tolerance: substituting a growing r into m^(4r+1) can be expensive. If the perturbation consumes the reference gap, rho>=1 and this guarantee does not apply.

For a certified subspace bracket beta<=Er<=v, the theorem also gives a ground-energy interval [beta-W rho^(2r+2),v]. This requires verifying the complete enrichment chain; it must not be applied to the current targeted 20- or 42-dimensional spaces as if they contained the full chain.

## Finite-fixture evaluation

For the five-pair Hamiltonians, B=binom(5,2)+5/5+eta=11+eta is a rigorous upper operator bound: same-side repulsion is at most binom(5,2), and the five reference hopping terms have total norm at most one. Use the independently certified b,u and the existing exact reference complement bound c. The exact rational calculation in `results/marginal_implicit_certificate/enrichment_order_bound.json` records rho and the first sufficient order for a target error of 1e-7.

At strength 1/100, rho is approximately 0.47917 for the cycle and 0.50786 for the mixed fixture. The first sufficient full-enrichment orders are 12 and 13, with guaranteed Ritz errors below 3.87e-8 and 4.54e-8 respectively.

The sufficient orders are conservative; on these small models a full enrichment can saturate the entire sector earlier. They do not describe the order needed by selective targeting. The new result is a controlled error statement within the declared perturbative region, not a general matter compiler.
