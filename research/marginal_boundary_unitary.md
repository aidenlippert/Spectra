# A certified contact-rotation limit and target-Hamiltonian transfer

The proposed covalent/ionic two-site rotation does not improve the existing
hopping-filter upper state. Its entire one-parameter family is now bounded
exactly, to about4.07×10^-10 per cut. The endpoint-RDM calculation also supports
fresh physical upper bounds for three target Hamiltonians with an interaction
absent from the original benchmark.

The best certified interval for the original U4,t1 Hubbard chain remains the
previous hopping-filter result: approximately[-0.611636,-0.558113972244] per
site at one million sites. The new rotation is a verified unsuccessful search
branch, and the transferred models currently have upper bounds only.

## Physical contact gate and local contraction

In the two-site occupation basis, define normalized singlets

\[
S=(|9\rangle-|6\rangle)/\sqrt2,\qquad
D=(|3\rangle+|12\rangle)/\sqrt2.
\]

For real x, put

\[
c=\frac{1-x^2}{1+x^2},\quad s=\frac{2x}{1+x^2},\qquad
V_x=I+(c-1)(SS^T+DD^T)+s(DS^T-SD^T).
\]

The verifier constructs all16 columns and checks exact orthogonality and both
spin-number conservation laws. Rational x gives a rational real unitary.

The source state is the existing physical H8 polynomial witness. Its left and
right two-site RDMs, each16×16, are recomputed by streaming the signed orbit
state. Trace, Hermiticity, and the independently recomputed source energy are
checked. No supplied density matrix is accepted as a physical certificate.

On the four sites surrounding a block cut, embed the gate on the middle two
sites. If ρ_R and ρ_L are the original endpoint RDMs, the energy increment is

\[
\Delta(x)=\operatorname{tr}\left[
(\rho_R\otimes\rho_L)(V_x^T H_4 V_x-H_4)\right].
\]

H_4 includes the onsite terms and all three nearest-neighbor bonds. The outer
onsite terms commute with the gate and cancel. This four-site neighborhood
contains every term changed by conjugation. The contraction covers all256
local occupation states without building a sixteen-site wavefunction.

For blocks of length8, gates on other cuts are disjoint from this changed
neighborhood. Their energy changes therefore add. Remote two-site RDMs are
also unchanged by a unitary acting on the complement. For q identical blocks,

\[
E_q=q e_8+(q-1)\Delta(x).
\]

This is a normalized physical circuit with no filter-normalization factor.
The statement is independent of any eigenstate or reflection assumption.

## Exact exclusion of the rotation family

The gate numerator is a quadratic matrix polynomial in x. Consequently

\[
\Delta(x)=\frac{a_0+a_1x+a_2x^2+a_3x^3+a_4x^4}{(1+x^2)^2}.
\]

The coefficients are recovered exactly from five rational evaluations, after
the degree bound has been established algebraically. For the actual H8 state
they are approximately

\[
(a_0,a_1,a_2,a_3,a_4)=
(0,-1.426445821887,4.135652325932,1.426445821887,1.282760682159).
\]

For the rational endpoint l=-0.110426865, the verifier constructs a3×3 Gram
matrix Q satisfying the exact polynomial identity

\[
\sum_{k=0}^4 a_kx^k-l(1+x^2)^2
=(1,x,x^2)Q(1,x,x^2)^T.
\]

Fraction-free integer PSD elimination accepts Q. This proves Δ(x)≥l for
every real x, including the rotation approached as x tends to infinity.
The rational parameter x=0.150424 gives

\[
-0.110426865\le\inf_x\Delta(x)
\le\Delta(0.150424)\approx-0.110426864593367.
\]

The exact bracket width is about4.0663×10^-10. A fresh contraction of the
existing linear filter gives its achieved shift-0.229107671571803, strictly
below this entire rotation family. Thus a better scalar angle cannot make the
proposed rotation competitive with that filter. The excluded family is this
single singlet-plane rotation on the fixed source state; other gates and other
building-block states remain open.

The Gram lower bound limits this variational family. It is not a ground-energy
lower bound for the unrestricted Hubbard problem.

## Transfer to a different local interaction

The source state and target Hamiltonian are now separated explicitly. The
source remains the validated H8 U4,t1 polynomial state. The target is

\[
H(U,t,V)=-t\sum_iT_i+U\sum_iD_i
+V\sum_i(n_i-1)(n_{i+1}-1).
\]

All parameters are bounded exact rationals, with U,t nonnegative. The current
transferred replay supports chain lengths divisible by8. This target family
includes a nearest-neighbor density interaction beyond the original onsite
Hubbard benchmark.

The code freshly computes the source expectations of total doublons and
internal density correlations. Since its original energy is e_8=4〈D〉−〈T〉,
the target block energy is exactly

\[
e_{\rm target}=U\langle D\rangle
-t(4\langle D\rangle-e_8)
+V\left\langle\sum_{i=1}^{7}(n_i-1)(n_{i+1}-1)\right\rangle.
\]

The unfiltered cross-cut density energy is the product
V〈n_R−1〉〈n_L−1〉. This vanishes for the accepted source, but is computed rather
than omitted from the formula. Cross-cut hopping has zero expectation before
the operation because both factors have fixed particle number.

Both the unitary and the linear filter are evaluated using the actual four-site
target Hamiltonian. For F=I−ηh_0, the local change is

\[
\operatorname{tr}\left[(\rho_R\otimes\rho_L)
\left(\frac{F^T H_4(U,t,V)F}{1+\eta^2}-H_4(U,t,V)\right)\right].
\]

Exact half occupation of every contact spin is checked before admitting the
filter branch. Its previously proved endpoint-preserving channel identity
then applies to entire remote even two-site RDMs, so the local target-energy
increments add at every merge. The code restores the uncentered physical
Hamiltonian before contraction: commutation with number alone would not justify
discarding a centered shift after nonunitary normalization.

Fresh million-site variational uppers, approximately, are:

| Target(U,t,V) | Fixed rotation | Fixed linear filter |
|---|---:|---:|
| (4,1,1/2) | -0.585425933742 | -0.600844010323 |
| (4,1,-1/2) | -0.501132046364 | -0.515383934166 |
| (3,2/3,1/2) | -0.372969825857 | -0.382603822649 |

The rotation uses x=0.150424 and the filter uses η=57277/250000. These are
legal physical upper states for each specified target, with all target energies
recomputed. Neither the state nor these parameters have been optimized for the
new targets. The table contains approximations, not outward-rounded standalone
decimal certificates; exact rational endpoints are in the replay artifact.
Each target and operation is also replayed at24 sites. No new-target lower
bound, exact ground energy, or generic molecular transfer is claimed.

## Independent validation and next attack

Direct sparse tests apply the contact operations to two and three entangled
four-site building blocks. They check the full Hubbard energy and normalization
against the endpoint contraction, and verify unchanged remote two-site RDMs.
The transfer test separately constructs the full12-site target energy for
U=3,t=2/3,V=1/2, for both operations, including a nonzero density-interaction
contribution. Thus the new interaction is exercised by an independent finite
state calculation.

```sh
python -S results/marginal_graded_hubbard8/discovery/boundary_unitary.py
```

The driver records exact family bounds,12 transferred chain replays, source
hashes, and the fixed source recipe under
`results/marginal_graded_hubbard8/boundary_unitary/`.

Final validation passed all444 regression tests in303.886seconds and all19
focused tests in15.201seconds. The standard-library replay completed in
25.347seconds, and all11 recorded source/input hashes matched afterward.

The next block-state objective can be derived exactly. At fixed η, define
w=η²/(1+η²) and Q_i=D_i−n_i/2+1/2. For a particle-hole symmetric source state,

\[
K_\eta=H_8+w\left[\frac12(T_L+T_R)-8(Q_L+Q_R)\right]
\]

satisfies

\[
e_8+g_\delta(\eta)=\langle K_\eta\rangle
+4w-\frac{2\eta}{1+\eta^2}.
\]

K_η preserves spin exchange and the existing particle-hole symmetry convention.
This converts optimization of the filtered block state into a precise modified
local Hamiltonian problem. No improvement from solving that new objective has
yet been demonstrated. Sharp accuracy at controlled cost, automatic lower
constraint discovery, and transfer to general chemistry remain open.
