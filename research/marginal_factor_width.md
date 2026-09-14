# Exact obstruction to three-configuration positivity pieces

The H6 reference gap cannot be compressed into a sum of positive semidefinite pieces each supported on at most three determinant configurations. An exact96-configuration dual witness proves this for the current reference and retained P. The obstruction is strong enough to defeat the entire scalar norm-shift transfer to hopping strength1/50, even when the perturbation norm is computed optimally on Q.

An explicit positive rank-one atom supported on four configurations rejects that dual witness. This supplies a concrete next constraint, but does not establish that four-coordinate pieces suffice. Configuration support here is a matrix-coordinate notion: it is not fermion body degree, number of interacting electrons, or orbital locality.

## The cone being tested

For a matrix A, the factor-width-k cone consists of sums

\[
A=\sum_j E_{S_j}^{T}A_jE_{S_j},\qquad
A_j\succeq0,\quad |S_j|\le k.
\]

The supports may overlap and may be any coordinate subsets; there is no restriction to a selected graph neighborhood or triangle dictionary. A dual matrix B is nonnegative on this cone if every principal submatrix of B of order at most k is PSD. These are established factor-width notions; see [the factor-width-rank paper](https://arxiv.org/abs/2405.11556) and [the study of k-PSD closures](https://arxiv.org/abs/2405.01208).

Here A=QH0Q-gamma Q in the determinant coordinates of the exact Sz=0 reference. The existing direct-spin factor certificate proves physical positivity at the target gamma. The new result excludes one restricted representation of that positive matrix, not the matrix's physical positivity.

## A compact dual witness with an exact global guarantee

Choose positive weights w on a subset S of Q. Define B to vanish outside S and set

\[
B_{ii}=w_i^2,\qquad
B_{ij}=-\operatorname{sgn}(H_{0,ij})\frac{w_iw_j}{k-1}\quad(i\ne j).
\]

On any subset of size r<=k, diagonal congruence by the weights leaves a unit-diagonal matrix whose off-diagonal absolute row sums are at most(r-1)/(k-1)<=1. That matrix is symmetric diagonally dominant with nonnegative diagonal, hence PSD. Zero weights outside S simply add zero rows and columns. Thus B is in the full dual cone, without enumerating all k-coordinate supports.

The exact normalized pairing is

\[
\frac{\operatorname{Tr}[(QH_0Q-\gamma Q)B]}{\operatorname{Tr}B}
=\frac{\sum_i(H_{0,ii}-\gamma)w_i^2
-\frac1{k-1}\sum_{i\ne j}|H_{0,ij}|w_iw_j}
{\sum_iw_i^2}.
\]

For the exported k=3 witness, replay gives **-0.030265072110328486 Ha** at gamma=-6.26489910104 Ha. Therefore A is outside the entire factor-width-three cone. The proof uses96 Hamiltonian source actions. Discovery used the explicit200/168 reference blocks and a numerical comparison eigenvector; the final checker uses only the exported integer weights and exact CAR actions. A negative result from this particular dual search would not prove cone membership.

Because the pairing is affine in gamma, the same witness proves an upper ceiling on any reference threshold that could admit such a representation:

\[
\gamma\le c
=-\frac{1519853240531305891726329722320283}
{241431867180473579748000000000000}
\approx-6.2951641731503285\;\mathrm{Ha}.
\]

This is stronger than a failed optimization or a failed chosen cluster cover. No number of three-coordinate PSD pieces can certify a higher threshold for this fixed H0, P, and basis.

## The Q-compressed norm escape also closes

For the strength1/50 connected hopping, exact CAR squares give ||DeltaH||<=1/25. The earlier saturating norm witness intersected P, so it did not establish optimality of ||Q DeltaH Q||.

The new physical witness has amplitudes[1,-1,1,1] on determinant states[627,630,633,636]. It lies entirely in Q and Sz=0, and its exact Rayleigh value under DeltaH is1/25. Consequently

\[
\|Q\Delta H Q\|=\|\Delta H\|=1/25.
\]

Only four source actions are required. The checker recomputes the square identities, physical support, exclusion from P, and saturation. It rejects the earlier witness[243,246,249,252] for this Q because two of those states lie in P.

Every factor-width-three reference proof followed by any valid scalar Q-norm shift therefore has threshold at most

\[
c-1/25\approx-6.3351641731503285\;\mathrm{Ha}.
\]

The independently replayed current Temple certificate places the Hs ground energy at least **0.001926057952290026 Ha above this ceiling**. Thus the shifted ceiling cannot exceed any valid upper witness, even an exact ground eigenvector. Improving the upper solver cannot repair this particular reference-gap pipeline.

The combined checker binds the same H0, P, sector, actual Hs-H0 perturbation, Q-norm saturation, and current Hs lower proof. Its replay includes the400-configuration current interval; only the standalone dual obstruction and standalone Q-norm proof have96- and4-source replay respectively. These counts must not be conflated.

## The next positive constraint is explicit

Within the96-state dual support, choose states[111,126,219,231] with signs[1,1,-1,1]. In coordinates rescaled by w, the dual quadratic form on this direction is exactly-2. Clearing the weight denominators gives an integer vector z on those four states.

The rank-one matrix zz^T is PSD, but its pairing with B is negative. After normalizing both matrices to unit trace, the exact pairing is approximately **-0.013403993686108635**. Hence the missing constraint

\[
z^TBz\ge0
\]

excludes this dual witness at factor width four. The new atom check uses four source actions; replay of its parent obstruction uses96. This is a verified separator and a candidate direction for a wider-cone search, not a complete factor-width-four decomposition of the reference.

## Two simpler compression attempts were weak

Deterministic minimum-degree symbolic elimination on the two exact reference graphs produced the following counts:

| Q block | Initial edges | Largest front, including pivot | Filled strict-lower entries | Dense strict-lower entries |
|---:|---:|---:|---:|---:|
|200|6,012|176|17,479|19,900|
|168|4,209|149|12,342|14,028|

This ordering saves about12% of strict-lower storage. It is not an optimal-ordering result or a treewidth lower bound. It gives no strong configuration compression on this fixture.

A separate numerical probe truncated the inherited LDL columns and bounded the remaining error by diagonal dominance. In the200-state block, retaining at most128 entries per column yielded only a candidate threshold near-6.34483 Ha, too low even before the perturbation shift. At160 entries per column the candidate improved to-6.26757 Ha, but retained18,851 of the original18,914 nonzero column coefficients. This is negligible coefficient compression. These are floating diagnostics, not accepted positivity certificates, and do not rule out other sparse factors.

## Replay and next boundary

```sh
python -S -m experiments.marginal_factor_width \
  --verify results/marginal_h6/factor_width3_obstruction/certificate.json
python -S -m experiments.marginal_operator_norm \
  --verify-q results/marginal_h6/hopping_q_norm_squares/certificate.json
python -S -m experiments.marginal_factor_width \
  --verify-transfer results/marginal_h6/factor_width3_q_transfer_obstruction/certificate.json
python -S -m experiments.marginal_factor_width \
  --verify-separator results/marginal_h6/factor_width4_separator/certificate.json
```

The implemented tests independently assemble B and its trace pairing with mixed-sign, non-unit coefficients, check all small principal inertias in a synthetic example, reject invalid supports and sectors, bind the actual perturbation, and distinguish a Q-supported norm witness from one overlapping P.

The live alternatives are wider factors, a different reference or retained space, coordinate transformations, and directional perturbation bounds that preserve more information than a scalar norm. The representation-size problem and general physical-marginal boundary remain open. This result identifies a genuine insufficient cone and an explicit constraint beyond it; it does not convert a finite obstruction into a universal lower bound on chemistry algorithms.

The full marginal regression suite passed **223 tests in225.387 seconds**. All five new standalone replays passed under `python -S` and match their saved receipts. The cumulative energy/witness ledger remains114; these new artifacts are tracked as obstructions, a norm certificate, and a separator.

A subsequent [constructive clique experiment](marginal_clique_gap.md) shows that support size alone does not settle the problem: all equal-amplitude support-four squares are insufficient for the current transfer, while a positive diagonal metric removes the corresponding necessary row obstruction. The general unequal-amplitude factor-width-four problem remains open.
