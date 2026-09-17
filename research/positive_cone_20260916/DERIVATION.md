# Exact positive-cone certificate for the eight-site Hubbard target

This derivation closes the **finite-model energy-width target** of `0.001t`.
It uses the full amplitude matrix, so it does not close the separate goal of
constructing and verifying a compact, scalable many-body representation.
All energies below are divided by the positive hopping amplitude `t`.

## 1. Exact physical model and sector

Label the open two-leg, four-rung ladder by `i=2r+s`, with
`r=0,1,2,3` and `s=0,1`. Its ten edges are

\[
\mathcal E=\{(0,1),(2,3),(4,5),(6,7),
(0,2),(1,3),(2,4),(3,5),(4,6),(5,7)\}.
\]

The dimensionless Hamiltonian is

\[
h=H/t=-\sum_{(i,j)\in\mathcal E,\sigma}
(c_{i\sigma}^{\dagger}c_{j\sigma}+c_{j\sigma}^{\dagger}c_{i\sigma})
+8\sum_{i=0}^{7}n_{i\uparrow}n_{i\downarrow}.
\]

We first work at `N_up=N_down=4`. The target is the **total energy width**,
not energy per site. This is a Hubbard lattice model; `t` is not a hartree
unit and the calculation is not molecular hydrogen H8.

The Hamiltonian commutes with total spin. At fixed total particle number
`N=8`, every allowed spin representation has integer spin and contains an
`M_S=0` state of the same energy. Therefore the minimum in the balanced
sector equals the minimum over all spin states at fixed `N=8`. This
argument does not compare energies at other particle numbers.

## 2. Particle-hole transformation and the amplitude map

Choose the bipartite signs

\[
\eta_{2r+s}=(-1)^{r+s}.
\]

Every hopping edge obeys `eta_i eta_j=-1`. Perform particle-hole conjugation
on down spin, with

\[
c_{i\downarrow}^{\dagger}\longmapsto\eta_i d_{i\downarrow},\qquad
n_{i\downarrow}\longmapsto1-n^h_{i\downarrow}.
\]

For distinct sites, anticommutation contributes a second minus sign to
the transformed hopping, so its matrix is unchanged. The interaction becomes

\[
8\sum_i n_{i\uparrow}(1-n^h_{i\downarrow})
=32-8\sum_i n_{i\uparrow}n^h_{i\downarrow}.
\]

There are four up particles and four down holes. Let

\[
d=\binom84=70,
\]

let `K` be the real symmetric one-spin hopping matrix on the four-particle
occupation basis, and let `N_i` be the diagonal one-spin occupation
projector for site `i`. A transformed state has an amplitude matrix
`C in C^(70 x 70)`. The exact Hamiltonian action is

\[
\boxed{\mathcal L(C)=KC+CK-8\sum_i N_i C N_i+32C.}
\]

In a general vectorization convention the right hopping term is `CK^T`;
here `K^T=K`. The Hilbert-space norm is the Hilbert-Schmidt norm,
`||C||_HS^2=Tr(C^dagger C)`. Consequently

\[
\frac{\langle\psi_C,h\psi_C\rangle}{\langle\psi_C,\psi_C\rangle}
=\frac{\langle C,\mathcal L(C)\rangle_{\rm HS}}
{\langle C,C\rangle_{\rm HS}}.
\]

The signed transformation was checked against independently constructed
Fock-space transitions on all 4,900 balanced configurations, and separately
on all 36 configurations of the four-site regression model. That enumeration
is a validation cost. The accepting checker reconstructs only the one-spin
basis and this amplitude map.

## 3. Why a ground eigenmatrix can be chosen positive semidefinite

The map `L` is self-adjoint under the Hilbert-Schmidt inner product and
preserves Hermitian matrices. For any complex matrix `C=X+iY`, where
`X` and `Y` are Hermitian,

\[
\|C\|_{\rm HS}^2=\|X\|_{\rm HS}^2+\|Y\|_{\rm HS}^2,
\qquad
\langle C,\mathcal L(C)\rangle
=\langle X,\mathcal L(X)\rangle+\langle Y,\mathcal L(Y)\rangle.
\]

The Rayleigh quotient is therefore a weighted average of the two Hermitian
Rayleigh quotients. A lowest eigenmatrix can be chosen Hermitian.

For Hermitian `X`, its energy numerator is

\[
q(X)=2\operatorname{Tr}(KX^2)
-8\sum_i\operatorname{Tr}(XN_iXN_i)+32\operatorname{Tr}(X^2).
\]

Replace `X` by `|X|=(X^2)^(1/2)`. Its norm, hopping contribution and constant
contribution stay unchanged. In an eigenbasis of `X`, with eigenvalues
`lambda_a`,

\[
\operatorname{Tr}(XN_iXN_i)
=\sum_{a,b}\lambda_a\lambda_b |(N_i)_{ab}|^2
\le\sum_{a,b}|\lambda_a\lambda_b| |(N_i)_{ab}|^2
=\operatorname{Tr}(|X|N_i|X|N_i).
\]

The coefficient of this expression is negative, so `q(|X|)<=q(X)`.
If `X` minimizes the Rayleigh quotient, `|X|` also minimizes it. Hence,
writing the ground energy of `h` as `e_0`, there is a matrix

\[
W\succeq0,\qquad W\ne0,\qquad\mathcal L(W)=e_0 W.
\]

This is the spin-reflection positivity structure associated with the
Hubbard model. It is established mathematics, rather than a new universal
positivity principle; see [Lieb, Two theorems on the Hubbard model (1989)](https://doi.org/10.1103/PhysRevLett.62.1201) and
[Boretsky, Cohn and Freericks (2017)](https://arxiv.org/abs/1712.02694).
The elementary argument above supplies exactly the property required here;
it needs neither uniqueness nor a quantitative spectral gap.

## 4. The accepting lower and upper bounds

**Certificate theorem.** If `C` is Hermitian and strictly positive definite
and the rational scalar `ell` satisfies

\[
C\succ0,\qquad D=\mathcal L(C)-\ell C\succeq0,
\]

then `ell<=e_0`. Indeed, products of two PSD matrices have nonnegative trace,
and self-adjointness of `L` gives

\[
0\le\operatorname{Tr}(WD)
=\operatorname{Tr}(\mathcal L(W)C)-\ell\operatorname{Tr}(WC)
=(e_0-\ell)\operatorname{Tr}(WC).
\]

Since `C` is strictly positive and `W` is nonzero PSD,
`Tr(WC)>0`. Division proves the lower bound.

The same matrix supplies a variational upper bound,

\[
\boxed{
\ell\le e_0\le u_C:=
\frac{\operatorname{Tr}(C\mathcal L(C))}{\operatorname{Tr}(C^2)}.
}
\]

These are **70-dimensional matrix positivity tests**. The theorem connects
them to the 4,900-dimensional physical Hilbert space. Positivity of an
arbitrarily reshaped residual would not be sufficient without Section 3.

Strict positivity cannot simply be relaxed to `C>=0`. For example, set
`K=diag(-1,1)`, `L(C)=KC+CK`, and `C=diag(0,1)`. Then `L(C)-2C=0`, although
the lowest eigenvalue of `L` is `-2`. The missing strictness allows the
candidate to be orthogonal to the ground PSD matrix.

For a general self-adjoint map there need not be any PSD ground eigenmatrix.
For example `L(X)=Tr(X)I-X` on two-dimensional matrices has `L(I)=I`, but
every traceless eigenmatrix has eigenvalue `-1`. This explains why the
Hubbard-specific sign structure is an essential hypothesis.

## 5. Proposal, rationalization and exact replay

Discovery starts from the identity matrix and applies the exact model's
amplitude map numerically in a Krylov eigensolve. It uses 4,900 amplitude
coordinates, with no previously computed state or energy endpoint. A
numerical solve is only a proposal.

For a positive candidate, a useful proposal for the lower endpoint is

\[
\lambda_{\min}\left(C^{-1/2}\mathcal L(C)C^{-1/2}\right),
\]

because `L(C)-ell C>=0` is equivalent to this eigenvalue being at least
`ell`. The code rounds the amplitude matrix to integers at scale `10^10`
and rounds the proposed lower downward. It then discards numerical
eigenvalue claims as evidence and runs exact checks.

An overall positive scale of `C` changes neither endpoint. With integer
`C` and `ell=p/q`, `q>0`, the checker constructs the integer matrix

\[
q\mathcal L(C)-pC
\]

directly from the ten hopping edges, `U/t=8` and the occupation labels.
Fraction-free elimination proves strict positivity of `C` and PSD of this
integer residual, including symmetry, pivot-sign, null-pivot and exact
division checks. The Rayleigh quotient is recomputed as a rational number.
The stored numerical energy and optimizer success status play no role in
acceptance.

If `mu=u_C`, a sufficient analytic version is

\[
\ell=\mu-\frac{\|\mathcal L(C)-\mu C\|_{\rm op}}
{\lambda_{\min}(C)}.
\]

The implementation uses the direct PSD gate instead of this potentially
looser norm estimate. The formula nevertheless exposes a possible scaling
problem: a small minimum eigenvalue of `C` can amplify residual errors.

## 6. Accepted numerical target

The cold candidate was accepted on its first rationalization attempt.
The exact endpoints are

\[
\ell=-\frac{30259232831}{10^{10}}=-3.0259232831,
\qquad
u_C=-\frac{151296140284988781488}{50000000000150236249}.
\]

Their exact difference is

\[
\boxed{
u_C-\ell=
\frac{238704658218758147090919}{500000000001502362490000000000}
\approx4.774093164360818\times10^{-7}<10^{-3}.
}
\]

Both integer positivity checks have rank 70. Fresh replay with Python's
standard library alone reproduced the exact endpoints and rejected all
five deliberately corrupted certificates. Five focused regression tests
passed, including the independent signed physical-Hamiltonian comparison.

## 7. What this settles and what it leaves open

The finite, fully coupled eight-site accuracy target is met. This proof
uses all excitation information present in the amplitude map instead of
replacing each patch by a single excitation threshold. It does not change
or invalidate the preceding exact limitation of the local proof family;
it uses a different family with more retained information.

The one-spin occupation basis has only 70 labels here, but its square
matrix contains all 4,900 amplitudes. For an even number `s` of sites,

\[
d=\binom{s}{s/2},\qquad d^2\sim\frac{4^s}{\pi s/2}.
\]

The current constructor, storage and exact PSD check therefore remain
exponential. No claim is made of enumeration-free discovery, tensor
compression, molecular transfer, a competitive general solver, dynamics,
inverse design or a resolution of the many-body problem.

The remaining mathematical requirement is a representation in which both
`C>0` and `L(C)-ell C>=0` can be constructed and checked without expanding
this full matrix, while keeping the energy interval useful. A compact
description of `C` alone would not supply that second inequality.
