# Exact extensive correlations from hopping filters

The physical upper state now correlates every neighboring pair of eight-site
blocks while retaining an exact energy formula independent of the number of
blocks. At one million sites, the verified interval for the half-filled uniform
open Hubbard chain at U=4,t=1 is

\[
-0.611636\le E_0/N\le-0.5581139722,
\]

with the displayed upper rounded outward. The exact width is approximately
0.0535220277555364 per site, compared with0.0821602575943402 for the preceding
uncorrelated eight-site tiling. This improves that width by about34.9%.

The construction is a restricted physical trial-state family. It supplies an
extensive upper bound, not a general representability oracle. Its best filter
parameter is also bracketed exactly, so further tuning within this fixed family
cannot close the remaining gap.

## One merge, without an eigenstate assumption

Let A and B be open clusters in normalized physical states with fixed total
spin populations. They can be internally entangled. Suppose each contact-site
spin occupation is exactly1/2. Define the hopping on the new bond by

\[
h=-\sum_\sigma(a_\sigma^\dagger b_\sigma+
b_\sigma^\dagger a_\sigma),\qquad
|\Psi_\eta\rangle=(I-\eta h)|\phi_A\rangle|\phi_B\rangle.
\]

The four directed spin-hopping outputs have distinct spin populations in A and
B, so their cross terms vanish. Each squared amplitude contributes1/4, giving
\(\langle h^2\rangle=1\). Every odd power of h flips the particle-number parity
of A. Consequently both \(\langle h\rangle\) and \(\langle h^3\rangle\) vanish.
The same sector argument eliminates terms with one h and H_A+H_B. Therefore

\[
\|\Psi_\eta\|^2=1+\eta^2,
\qquad
E(\eta)=e_A+e_B+
\frac{-2\eta+\delta\eta^2}{1+\eta^2}.
\]

To compute δ, the verifier checks this local operator identity by exact CAR
normal ordering, for each endpoint spin:

\[
cHc^\dagger+c^\dagger Hc
=H+(c_j^\dagger c+c^\dagger c_j)
  +U(n_{\bar\sigma}-2D).
\]

Here j is the adjacent site inside the same cluster. This identity holds as an
operator statement. It never replaces an expectation of nH by a product of
expectations. In particular, the trial state need not be an eigenstate.

Writing T_A,T_B for the positive spin-summed internal hopping next to the
two contacts, and D_A,D_B for contact doublon expectations, it follows that

\[
\delta=\frac{T_A+T_B}{2}+U(1-2D_A-2D_B).
\]

## Why every cut can be correlated

Let O_A be an even observable supported away from A's active contact site.
For the endpoint observables used here, it also preserves spin populations.
The partner's half occupations and the orthogonality of its directed hopping
outputs give

\[
\langle hO_Ah\rangle
=\frac12\sum_\sigma
 \langle a_\sigma O_Aa_\sigma^\dagger+
 a_\sigma^\dagger O_Aa_\sigma\rangle
=\langle O_A\rangle.
\]

The last equality uses disjoint support and the CAR anticommutator. The linear
filter terms vanish by particle-number sectors. Dividing by1+η² shows that
O_A's normalized expectation is unchanged. The same holds for B.

Thus a merge leaves the remote endpoint occupations, doublons, and internal
endpoint hopping unchanged. For eight-site building blocks, each outer
two-site region is disjoint from the joining contact. These are exactly the
observables needed at the next merge.

After joining A to B, the new factor is the entire cluster AB. It still has
fixed **total** spin populations. The induction does not require the old
constituent particle numbers N_A and N_B to remain fixed. At the next step the
product factors are AB and C, and the directed output sectors are classified
using the total charges of those two factors. This establishes the induction
even though the merged cluster is internally entangled.

For q original blocks and a common η, define

\[
g(\eta)=\frac{-2\eta+\delta\eta^2}{1+\eta^2}.
\]

Then the exact correlated trial energy and norm are

\[
E_q=q e_8+(q-1)g(\eta),\qquad
\|\Psi_q\|^2=(1+\eta^2)^{q-1},
\]

where the input blocks are normalized. The filters act on disjoint contact
pairs, and are even operators, so their product can equivalently be specified
without choosing an order. The proof uses sequential merging.

The implementation stores the normalization as a positive base and an integer
exponent. It does not construct its exponentially large numerator or a global
wavefunction. An unfiltered2-,4-, or6-site physical remainder can be appended;
the last hopping expectation then vanishes by fixed particle number of the
correlated cluster and remainder.

## Verified H8 data and the filter-family limit

The original physical H8 polynomial witness has

\[
e_8\approx-4.235805939245279.
\]

Its exact endpoint occupations are1/2 for each spin. The computed endpoint
values are approximately

\[
D_L=0.071694272264724,\quad D_R=0.071694272263592,
\]

\[
T_L=1.282760682161060,\quad T_R=1.282760682156884,
\]

giving δ≈4.135652325932441. The tiny left/right differences are retained;
reflection symmetry of the trial coefficients is not assumed.

The accepted rational parameter is

\[
\eta=57277/250000=0.229108,
\quad g(\eta)\approx-0.229107671571803.
\]

The shift is the Rayleigh quotient of

\[
M=\begin{pmatrix}0&-1\\-1&\delta\end{pmatrix}
\]

on the vector(1,η). An exact2×2 PSD check proves

\[
-0.229107672\le\inf_{\eta\in\mathbb R}g(\eta)
\le g(57277/250000).
\]

The exact bracket width is approximately4.282×10^-10 per cut. Allowing different
real η_j at different cuts does not escape the bracket, because δ is inherited
at every merge and the shifts add. This is a limit of the fixed physical
upper-state family, not a lower bound on the unrestricted ground energy.

The fresh chain replays give:

| Sites | Correlated upper energy/site, approximately |
|---:|---:|
| 8 | -0.529475742406 |
| 16 | -0.543794971879 |
| 18, including a two-site remainder | -0.529397037489 |
| 24 | -0.548568048370 |
| 64 | -0.554534393984 |
| 1,000,000 | -0.558113972244 |

The asymptotic density of this trial family is (e_8+g)/8. The finite open-chain
formula retains the missing final cut; no thermodynamic extrapolation is used
to certify a finite chain.

## Verification and remaining frontier

The code first independently replays the H8 moment upper. It then constructs
the same polynomial state by Horner recurrence in the existing signed orbit
basis, checks its energy against that moment quotient, and streams orbit
members to compute endpoint observables. It stores1,239 orbit amplitudes and
retains the existing sparse-action limits. It constructs neither a charged
cluster wavefunction nor a sixteen-site sector.

An independent direct CAR test uses a non-eigenstate four-site singlet as its
building block. Two- and three-block states are filtered on every cut, using
only72 and864 sparse amplitudes respectively. The tests verify the exact norm,
the full Hubbard energy including onsite terms, fixed total spin populations,
and unchanged remote endpoint observables. This checks the induction after a
cluster has become entangled, rather than only testing isolated pairs.

The exact driver is

```sh
python -S results/marginal_graded_hubbard8/discovery/hopping_filter.py
```

It recomputes endpoint data and the six-site lower proof, brackets the filter
family, and checks all reported finite chains. Inputs, receipts, and source
hashes are in `results/marginal_graded_hubbard8/hopping_filter/`.

All437 regression tests passed in293.211seconds. The initial five focused tests
passed in1.785seconds; the full run includes their final names and contents.
The fresh standard-library certificate replay completed in59.224seconds,
including the six-site lower PSD proof and six independently reconstructed
upper states. All11 recorded source/input hashes matched afterward.

The remaining accuracy gap now requires changing the building-block state,
the boundary operation, or the lower-certificate family. Further scalar filter
tuning has a proved limit. General molecular transfer and accuracy-versus-cost
guarantees remain open. The filters define valid normalized physical states;
this does not establish an efficient postselection or laboratory preparation
procedure.

Follow-up: [the contact-rotation audit](marginal_boundary_unitary.md) now
certifies that the proposed rational singlet/ionic rotation cannot beat this
filter on the fixed source state. The endpoint-RDM contraction also gives
physical upper-state transfer to specified nearest-neighbor density-interaction
targets. The next derived objective optimizes the source block for the energy
of the complete filtered chain.
