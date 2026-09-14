# Joint polynomial charge metrics

The newer [joint-metric and vanishing-metric constructions](marginal_joint_metric.md) remove every singleton positive atom on both H6 Hamiltonians. Exact gap and energy replays pass at −6.24 for a jointly optimized degree-four metric and −6.264 for a degree-six metric that vanishes on the valence space. The latter removes the global projector from compilation. The existing energy widths and 17 response directions are preserved. **353 full marginal tests pass.** Number-ideal lifting, discovery, and energy-response work remain finite-sector computations. The earlier fixed-metric obstructions below remain valid in their stated scope.

The joint transition dependence can be retained in one exact occupation polynomial. On the finite localized H6 model, a degree-four charge metric certifies the complement threshold **−6.24 Ha**, improving the preceding −6.264 Ha certificate. This clears the independent-ratio-envelope obstruction described in `marginal_charge_tail.md`.

The first accepted proof still contains 159 singleton occupation atoms. It is a new exact algebraic route, not evidence that many-body complexity has disappeared.

## The representation and inequality

Write the spatial charge as q_i=n_(i,alpha)+n_(i,beta)−1. Let v(q) be a rational polynomial, strictly positive on the ionic complement Q of the complete singly occupied valence space P. For the fixed-spin sector, compile

\[
K_\gamma(n)=(V(n)-\gamma)v(q)
 -\sum_g I_g(n)\sigma_g A_g(n)v(q+\Delta_g)Q(q+\Delta_g).
\]

Here I_g is the source occupation event, A_g is the coherently grouped CAR amplitude, and sigma_g is its checked fixed sign on that event. The exact target projector is

\[
Q(q')=1-\prod_i(1-(q'_i)^2).
\]

This expression preserves the dependence between all transition amplitudes, source events, target charges, and metric values. It does not replace each metric ratio with a separately attained maximum. No division by v is needed during compilation.

If v>0 and K_gamma≥0 on Q, every weighted Q row is diagonally dominant at gamma. Diagonal similarity followed by Gershgorin bounds the Hermitian complement spectrum below by gamma. This is a sufficient positivity certificate for this Hamiltonian and retained space; it is not a description of the full physical marginal cone.

The metric may be nonzero on P. The source inequality is asserted only on Q, while target-P transitions are explicitly removed. Negative values of K on P are therefore permitted and tested.

## Exact polynomial proof

Occupation multiplication uses n_i²=n_i. Any monomial containing more occupied alpha or beta modes than the fixed population is zero on the physical slice and is discarded exactly. The checker reconstructs

\[
p=b+\sum_a c_a I_a +(D-1)\sum_a d_a I_a
 +(N_\alpha-3)L_\alpha+(N_\beta-3)L_\beta+\epsilon,
\]

for p=v and p=K_gamma. Indicator coefficients c_a,d_a are nonnegative rationals. D counts doublons. The number multipliers have unrestricted signed rational coefficients. On Q, D≥1, so the certified lower bound is b−||epsilon||_1. This residual is recomputed in exact integer arithmetic from the supplied Hamiltonian and metric.

The energy verifier binds the recipe to the actual Hamiltonian, requested threshold, and complete valence retained space. Recipes containing hidden Hamiltonians or thresholds are refused. A changed Hamiltonian requires its own checked proof or a separately checked transfer theorem.

## Discovery and accepted finite proof

Numerical metric discovery used all 380 Q rows and reflection-tied charge features. With a constant term allowed, degree four and 48 independent features suffice. Requiring v(0)=0 forced a numerical search to degree six and 73 independent features, the entire reflection-symmetric charge-value space in this fixture. The latter is explicitly not compression.

The degree-four coefficients were rounded to denominator 10^8. Independent exact determinant actions agree with the compiled numerator on all 400 spin-sector occupations, including P sources. On the 380 Q sources, the metric minimum is 0.03775165 and the actual weighted-row minimum is −6.236835969910392 Ha.

Compilation uses 1,140 grouped transitions, 290,817 polynomial products, and 7,992,302 coefficient products, with 1,764 peak monomials. The final numerator has degree six after fixed-population pruning. The symbolic replay makes no determinant actions or occupation/charge endpoint loop. These statements concern replay, not the enumerative numerical discovery.

The first exact proof has:

| Quantity | Metric positivity | Numerator positivity |
|---|---:|---:|
| Positive indicators | 252 | 295 |
| Charge-localizer indicators | 56 | 103 |
| Number-multiplier terms | 1,290 | 1,364 |
| Singleton positive indicators | 0 | 159 |
| Certified lower bound | 0.009999999166 | 0.0009999988884280846 |

The numerator singleton atoms each fix three occupied alpha and three occupied beta orbitals. The fixed populations determine every other occupation. Their presence encodes individual determinants even though the verifier does not enumerate a list of rows. The JSON gap certificate is 790,440 bytes, compared with 487,364 bytes for the preceding charge-spin gap file; these sizes include different metric/proof layouts and are not an asymptotic comparison.

The first broad coefficient LP timed out. A reduced 23,007-column LP solved both positivity proposals. One numerator coefficient was −3.07e−12 numerically; an explicitly recorded proposal projection set it to zero. The exact residual checker, rather than LP status or tolerance, accepted the resulting rational certificate.

## Transfer and remaining pressure points

The corrected exact row diagnostic uses Q sources and Q targets in both Hamiltonians. The same metric has the same minimum −6.236835969910392 Ha under the checked 1/50 hopping perturbation. An earlier diagnostic included P sources and is invalidated; its apparent failure at state 1638 does not apply to Q. Reusing the original numerator decomposition under perturbation fails its conservative residual bound. A fresh exact decomposition now passes for the changed Hamiltonian at −6.24 Ha, with 271 positive indicators, 126 charge-localizer indicators and 1,364 number-multiplier terms. The metric and its positivity proof are unchanged.

Removing all degree-six positive atoms makes the current degree-four/localizer basis numerically infeasible. Adding every degree-five positive indicator does not restore feasibility. These are basis-dependent numerical outcomes, not exact impossibility theorems. Larger bases excluding singleton atoms have produced solver errors; an error is not evidence of infeasibility.

### Exact obstruction to removing singleton atoms from this cone

A corrected 400-coordinate LP, with the constant bound free and zero ideal columns removed, supplies an exact rational separating functional. Its nonzero support consists of twelve occupations of weight +1/11 and one occupation, state 819, of weight −1/11. The positive occupations have two doublons; the negative occupation has three.

Independent replay verifies normalization, all **86,509** admissible nonsingleton indicators of degree at most six, and all **9,909** nonempty charge-localizer indicators of degree at most four. Every allowed atom has nonnegative functional value. Number identities vanish because every support occupation has the required spin populations. Every absolute monomial moment is at most one, so the functional also controls the coefficient-L1 residual allowance.

Nevertheless,

\[
L(K_{-6.24})=-4487406304359711214887/55000000000000000000000
\approx-0.08158920553381294.
\]

Consequently this entire restricted cone cannot certify a nonnegative numerator lower bound for the current metric and threshold, even allowing arbitrary numbers of its atoms and the checker's bounded residual. The displayed number is an exact upper bound on the best certifiable constant; no matching exact primal optimum is claimed. It does not refute the successful proof that includes singleton atoms, a different metric, or a different positivity cone.

The witness identifies a concrete missing physical constraint:

\[
L((D-2)^2)=-1/11<0.
\]

The polynomial is a square and is nonnegative on every occupation. Adding its products with indicators of degree at most two stays within degree six. That enlarged finite LP improves the numerical best bound to about −0.06611061 but still falls below zero; one square family is insufficient. General polynomial-square separation is the next tested extension. The signed functional is a proof-family adversary, not a physical state or an electronic energy witness.

Twelve rounds of quadratic-moment negative-eigenvector cuts improve the numerical bound from −0.08158921 to −0.06785826, still negative. A separate family of triple-site number/doublon squares and their degree-two indicator localizers reaches −0.07553174. These runs supply directions for stronger positivity, not a successful singleton-free gap or a limit on all polynomial squares. The positive certificates above are unchanged.

The original energy proof accepts the new gap while retaining its 17 response directions and exact interval width 7.416980322924775e−12 Ha. The perturbed energy proof also accepts the direct new gap, retaining 17 response directions and width 3.9293895949815965e−12 Ha. Its nearby spin-symmetric width is 3.331192213308058e−14 Ha; the larger input-Hamiltonian interval includes the checked symmetry approximation allowance. These do not improve the intervals themselves: response accuracy and the approximation error already control them. No fresh upper-witness or response discovery is claimed; the existing objects were rechecked with the new gaps.

Still open: bounded-degree proof families under increasing system size, scalable upper-witness and response construction, molecular-family transfer beyond this fixture, general marginal representability, and the thermodynamic/kinetic/synthesis layers needed for materials design. Ground-state electronic energy alone does not solve those layers.

Machine-readable acceptance, independent replay paths, subsequent basis experiments, and current test results are recorded in `results/marginal_h6/polynomial_metric/progress.json`.

Validation: **345 marginal regression tests passed in 284.943 seconds**. Both new gap certificates, both integrated energy certificates and the restricted-cone separator pass independent standard-library-only replay. The gap and separator receipts match their exporter receipts exactly. Six new tests cover exact row reconstruction, target-valence exclusion, physical-slice reduction, positivity/sign refusals, energy/reference binding to the actual Hamiltonian, threshold and complete retained space, and the separator's normalization, atom signs and residual allowance.
