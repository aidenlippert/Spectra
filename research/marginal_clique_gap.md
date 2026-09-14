# Four-coordinate supports also need amplitude geometry

A constructive H6 clique-square certificate now replays exactly, but its complementary threshold is only **-7.553096043350333 Ha**, far below the transfer requirement. A separate one-row proof shows that adding more equal-amplitude four-coordinate squares cannot repair this family. Unequal coefficient ratios or a changed representation are necessary.

Positive diagonal rescaling removes that necessary row obstruction on both reference blocks. The rescaling is a concrete proposal for the next representation; it is not yet a complete gap certificate.

## The constructive cone

For a sign-consistent clique C of r=3 or4 determinant coordinates, choose signs sigma_i in{+1,-1} such that sign(H_ij)=sigma_i sigma_j on every clique edge. The atom

\[
\lambda_C\,\sigma_C\sigma_C^T,\qquad\lambda_C\ge0,
\]

is positive semidefinite. Its nonzero amplitudes have equal magnitude. The constructor allocates such atoms while keeping every edge budget nonnegative:

\[
\sum_{C\supset\{i,j\}}\lambda_C\le|H_{ij}|.
\]

The unused edge coefficient gives another positive square with support two. Relative to those edge squares alone, an r-coordinate clique saves(r-2)lambda_C from each participating diagonal. Thus the exact diagonal threshold is

\[
\gamma\le\min_i\left[H_{ii}-\sum_{j\ne i}|H_{ij}|
+\sum_{C\ni i}(|C|-2)\lambda_C\right].
\]

Replay reconstructs every block from the actual Hamiltonian, checks the clique signs and exact rational edge budgets, and checks this diagonal margin. No large LDL positivity factor is used. Complete configuration coverage still remains.

The LP only proposes nonnegative atom weights. Export rounds them to a common integer scale and repairs overspent edges by decreasing incident weights. Such decreases cannot overspend another edge. The final threshold is recomputed exactly; the numerical objective is never used as the accepted bound.

## Finite construction result and solver limits

Two all-column LP attempts reached their120-second limits without returning a primal proposal. Those failures do not establish infeasibility. The subsequent native HiGHS column-generation run retained its LP basis while adding batches of512 dual-priced cliques. It generated an exact certificate even though numerical optimization stopped before complete pricing convergence.

| Reference block | Candidate cliques | Selected LP cliques | Positive exported cliques | Exact block threshold (Ha) |
|---:|---:|---:|---:|---:|
|200|64,820|13,312|807|-7.553096043350333|
|168|32,013|10,554|560|-7.108785607496666|

Candidate generation used Hamiltonian edges of magnitude at least1e-8. Smaller edges still remain in the exact leftover-edge proof. No claim is made that this dictionary represents the entire factor-width-four cone, or that the returned LP point optimizes even this dictionary.

The certificate contains1,367 clique atoms and9,867 nonzero leftover edge squares, occupies458,421 bytes, and replays368 Q source actions while referencing400 spin-sector configurations. The recorded construction took186.51 seconds. The native run stopped after91.37 and92.16 solver seconds on the two blocks because its cumulative solver limit was being reduced twice; that timing setup was corrected afterward. No speedup or reproducibility of the time-limited LP trajectory is claimed. The saved exact certificate is independent of that discovery issue.

The gap is too low to enter the current positive-denominator transfer. It is recorded as a valid but insufficient complement bound, not as a new ground-energy interval.

## An exact obstruction to all equal-amplitude support-four squares

The obstruction covers more than the selected sign-consistent cliques. Suppose

\[
A=\sum_\alpha\lambda_\alpha z_\alpha z_\alpha^\dagger,
\quad\lambda_\alpha\ge0,
\]

where each z_alpha has at most k nonzero coordinates and those coordinates have equal magnitude. Arbitrary signs, phases, overlaps, and coefficient cancellations are allowed. For any row i, the triangle inequality gives

\[
\sum_{j\ne i}|A_{ij}|
\le\sum_{\alpha:i\in\operatorname{supp}z_\alpha}
\lambda_\alpha(k-1)|z_{\alpha,i}|^2
=(k-1)A_{ii}.
\]

Consequently, representing QH0Q-gamma Q this way requires

\[
\boxed{\gamma\le H_{0,ii}-\frac1{k-1}\sum_{j\in Q,\,j\ne i}|H_{0,ij}|.}
\]

At k=4, the exact row for determinant111 yields

\[
\gamma\le-6.644885509747556\;\mathrm{Ha}.
\]

Only one Hamiltonian source action is needed for this standalone obstruction. Combining it with the previously proved optimal compressed norm ||Q DeltaH Q||=0.04 places every scalar-shifted threshold at most-6.684885509747556 Ha. That is **0.3516473945495171 Ha below** the independently certified current Hs ground lower. Neither more flat four-coordinate atoms nor a more accurate upper witness can repair this fixed-reference pipeline.

This is an equal-amplitude restriction, not a general factor-width-four obstruction. A test makes that distinction explicit: the positive rank-one matrix[1,4][1,4]^T is already a two-coordinate square, yet violates the support-four equal-amplitude row condition at threshold zero. It requires unequal amplitudes.

## A diagonal metric removes the necessary obstruction

Let W=diag(w_i) have strictly positive entries. Flat squares in WH0W correspond to unequal-amplitude squares with coefficients proportional to1/w_i in the original coordinates. The transformed target is

\[
W(H_0-\gamma I)W=WH_0W-\gamma W^2.
\]

The identity coefficient therefore becomes the metric W²; treating it as gamma I would certify the wrong problem.

The transformed necessary row ceiling is

\[
c(W)=\min_i\left[H_{0,ii}-\frac{\sum_{j\ne i}|H_{0,ij}|w_j}{3w_i}\right].
\]

The numerical proposal uses the lowest eigenvector of diag(H0_ii-gamma0)-|offdiag H0|/3, takes positive absolute values, and rounds the normalized components to positive integers. Exact replay then recomputes c(W), without relying on the floating eigenvalue.

| Q block | Original uniform ceiling (Ha) | Rescaled ceiling (Ha) | Reference target (Ha) |
|---:|---:|---:|---:|
|200|-6.644885509747556|-6.127741542713512|-6.26489910104|
|168|-6.443756774655111|-5.951608754347920|-6.26489910104|

Both necessary row conditions now pass at the reference target. These numbers are upper ceilings permitted by a necessary condition, **not lower spectral bounds for H0**. The rescaling initially supplied only this necessary gate. The subsequent construction and exact family obstruction below now show that these particular metrics cannot produce the required balanced-clique threshold.

For the different auxiliary matrix diag(H0)-|offdiag H0|/3, these same values are rigorous weighted-Gershgorin lower bounds. Thus that auxiliary comparison matrix passes positivity at gamma0: the particular comparison-based dual family that obstructed factor width three cannot also supply a width-four obstruction here. More general width-four dual witnesses remain possible. The auxiliary matrix must not be substituted for the physical Hamiltonian.

This locates the next construction problem precisely: allocate wider positive directions in a learned diagonal metric, preserve the gamma W² target exactly, and test whether their amplitudes can close the gap with controlled proof size. A successful necessary gate does not replace that construction.

## Replays and remaining scope

```sh
python -S -m experiments.marginal_clique_gap \
  --verify results/marginal_h6/clique_width4_gap/certificate.json
python -S -m experiments.marginal_clique_gap \
  --verify-flat results/marginal_h6/flat_width4_obstruction/certificate.json
python -S -m experiments.marginal_factor_width \
  --verify-transfer results/marginal_h6/flat_width4_transfer_obstruction/certificate.json
python -S -m experiments.marginal_clique_gap \
  --verify-metric results/marginal_h6/clique_metric_probe/certificate.json
```

The last command verifies only the necessary metric condition. The transfer obstruction replays a complete existing400-configuration energy proof in addition to the one-row obstruction and Q-norm saturation. These proof scopes and source counts remain separate.

The broader representation-size and universal marginal-boundary problems are unresolved. This turn supplied a constructive but weak wider-cone certificate, a proof that flat amplitudes are insufficient here, and explicit rational metric weights that remove that particular obstruction.

The full marginal regression suite passed **230 tests in272.506 seconds**. All four new standard-library checks match their saved receipts. One is a necessary-gate diagnostic, not a positivity certificate. The cumulative accepted energy/witness ledger remains114.


## Actual weighted construction and an exact family ceiling

The constructor now accepts the rational diagonal metrics and verifies

\[
WH_0W-\gamma W^2
\]

as an exact sum of balanced clique squares, leftover edge squares, and nonnegative diagonal residuals. With transformed baseline \(b_i=A'_{ii}-\sum_{j\ne i}|A'_{ij}|\) and \(m_i=W_{ii}^2\), the exported threshold is \(\min_i[b_i+\sum_{C\ni i}(|C|-2)x_C]/m_i\). The metric coefficient is retained in both the numerical LP and exact replay. A nonconstant-metric rank-one test independently reconstructs the original matrix and rejects an overstated threshold.

The bounded construction exports a valid but weak common lower bound **−8.449457060735654 Ha**. It uses809 and856 clique atoms in the200/168 blocks, with6012 and4209 nonzero leftover edge squares. Both searches stop at the33-round pricing budget with negative unpriced reduced costs; their restricted LP statuses being “Optimal” do not mean the complete dictionary is optimized. Total recorded proposal/export time is275.644 seconds. This result alone is not evidence of cone insufficiency.

The independent dual resolves that ambiguity. For transformed edge capacities \(c_e=|A'_e|\), nonnegative node and edge weights \(y_i,z_e\) satisfying

\[
\sum_{e\subset C}z_e\ge (|C|-2)\sum_{i\in C}y_i
\]

for every balanced triangle and quadruple give the exact upper bound

\[
\gamma_{\rm packing}\le
\frac{\sum_i b_i y_i+\sum_e c_e z_e}{\sum_i m_i y_i}.
\]

Integer weights need no approximate normalization. Numerical proposals are rounded, every violated clique inequality is repaired by increasing an edge weight, and replay independently checks every inequality. The full nonzero-edge dictionary contains64,820 and32,013 cliques; it happens to equal the cutoff dictionary on these blocks.

The168-state dual proves

\[
\boxed{\gamma_{\rm packing}\le
-\frac{22695281579233455165226345131}
{3338288270638084000000000000}
=-6.798478663107022\;\mathrm{Ha}.}
\]

This is0.5335795620670224 Ha below the reference target. It excludes the entire balanced-clique packing family in these fixed metrics, regardless of further LP optimization. It does not exclude general factor-width-four decompositions, other metrics, or cancellation between differently signed squares.

The200-state numerical dual times out and needs a large exact feasibility repair, giving only−4.1928118133660055 Ha. That weak block does not undermine the independent168-state obstruction: a common threshold cannot exceed either block's valid ceiling. The proposal file's `solver_seconds` field includes its subsequent rational repair and verification; it is not a pure solver benchmark.

## The obstruction reduces to a local edge-cover proof

The decisive dual has one nonzero node weight at determinant3252 and41 nonzero edge weights, all equal after a common rescaling. Thus it is a combinatorial cover: every balanced triangle containing that row must contain at least one covered edge, and every balanced quadruple must contain at least two.

Standalone replay discovers the complete49-neighbor Q row neighborhood directly from the Hamiltonian. It checks302 triangles and410 quadruples. Cliques outside the active row have zero left side in the dual inequality and require no inspection. The same exact ceiling now uses **50 source actions and191 referenced determinants**, compared with368 source actions for the full dual replay. The certificate is247,925 bytes, mostly the shared Hamiltonian. No component list, complete sector enumeration, dense factor, or numerical optimization is needed in this local replay. This is a compressed obstruction, not a compressed positive energy certificate.

The local proof also applies to any positive extension of the declared metric weights outside those50 states. It does not allow changing the50 local weights. The coefficient vector was discovered using the explicit full-block LP, so discovery has not acquired this smaller frontier merely because verification has.

There is an explicit way out of the excluded family. Write its unnormalized dual matrix as

\[
B_{ii}=y_i,\qquad
B_{ij}=\tfrac12\operatorname{sign}(H_{ij})(z_{ij}-y_i-y_j).
\]

Six covered edges do not touch the active row. On each such pair, the diagonal entries of B vanish and its off-diagonal entry is \(\operatorname{sign}(H_{ij})/2\). The positive square with transformed amplitudes \((1,-\operatorname{sign}(H_{ij}))\) therefore has exact pairing−1 with B. These directions are recorded in `clique_metric_star_obstruction/cancellation_escape.json`. They show that permitting cancellation can invalidate this particular dual immediately. They do not prove that the enlarged family supplies the H6 gap; the prior general factor-width-three obstruction remains valid.

The next construction should allow independent local square directions and a diagonally dominant exact residual, including both signs of edge squares. A dense PSD residual would merely restore the old full-factor dependency. Whether a small such dictionary reaches the target remains open.

```sh
python -S -m experiments.marginal_clique_gap \
  --verify results/marginal_h6/clique_metric_width4_gap/certificate.json
python -S -m experiments.marginal_clique_dual \
  --verify results/marginal_h6/clique_metric_width4_dual/certificate.json
python -S -m experiments.marginal_clique_dual \
  --verify-star results/marginal_h6/clique_metric_star_obstruction/certificate.json
```

All three independent standard-library replays match their saved fields. These are a gap certificate and two family-obstruction certificates; the accepted energy/witness ledger remains114.

The final regression suite passed **236 tests in252.183 seconds**. The new focused clique/dual suite passed13 tests; independent replays use only the standard library. The preceding230-test record above describes the earlier unweighted step.

A subsequent numerical dual probe adds the constraints z_ij <= 2(y_i+y_j), corresponding to allowing both signs of residual edge squares. On the168-state block, HiGHS reports an optimum−6.74641136040593 Ha with maximum reported primal violation about1.71e−13. The ceiling is higher than the previous−6.798478663107022 bound because the primal square family is enlarged. This probe is **not certified**. Its168 node and4209 edge weights are preserved in `clique_metric_star_obstruction/both_sign_numerical_probe.json` for rational repair and exact checking. No accepted obstruction claim depends on this numerical result.

The previously noncertified both-sign proposal is now exact-checked in [the signed-clique report](marginal_signed_clique.md). That report also gives a verified square outside the balanced dictionary and two actual signed-atom complement certificates, still below the required threshold.
