# Joint occupation algebra: repaired local bounds, remaining global cost

The joint occupation compiler repairs the saved failing H6 branch exactly. A one-electron/one-hole simplex extension then certifies the entire ionic complement at **−6.264 Ha** with **335 tree nodes and 168 shared leaves**. Both original and perturbed energy intervals remain verified. This is a finite certificate improvement: the simplex leaves still evaluate **400 physical occupation endpoints**. General representability and scalable molecular compilation remain open.

## The specific dependence that was lost

The previous source-indicator bound treated the transition amplitude and positive charge-metric ratio with separate worst cases. The saved branch at mask1535, bits1504 leaves modes9 and11 unresolved with one beta electron. Both physical row inequalities pass, but the independent envelope gives −6.384536849053956 Ha and misses the −6.264 threshold.

Set n9=x and n11=1−x in the Boolean quotient algebra

\[
\mathbb Q[x]/(x^2-x).
\]

The compiler retains source event, coherent CAR amplitude, target-valence exclusion, and metric ratio together. The complete weighted row becomes R(x)=c0+c1*x, with

- c0≈−5.544518805707165 Ha;
- c1≈−0.3857306751986948 Ha;
- min R=c0+min(0,c1)≈**−5.93024948090586 Ha**.

The endpoints match independently computed physical rows exactly as rational numbers. Thus this branch obstruction is repaired, rather than merely bypassed by selecting a looser target. The saved receipt contains the exact fractions.

Absolute values and positive rational exponentials in this algebra use two scalar endpoint evaluations. Keeping their dependencies together is exact, but does not eliminate that local endpoint work.

## From a line to an occupation simplex

For k unresolved same-spin modes containing one electron, use the vertex algebra

\[
x_i x_j=\delta_{ij}x_i,\qquad \sum_i x_i=1,\qquad n_i=x_i.
\]

For one hole, n_i=1−x_i. Each element has k endpoint coordinates; addition, multiplication, absolute value, and metric exponentials act coordinatewise. This supports a k-point physical simplex instead of the full 2^k Boolean cube. It does not compress those k physical points further.

The compiler processes each transition group jointly over the coordinates and caches its target-projector/metric kernel by the vector of charge changes. It constructs neither determinant integer labels nor determinant actions. Source-valence coordinates receive B=sum|h_word|. This sentinel preserves the minimum over Q because R_Q(s)≤H_ss≤B and the chart must contain at least one Q state. Inactive source-event coordinates multiply away formal target-charge expressions outside an allowed transition.

## Construction and replay have separate costs

A new optional construction strategy splits an ineligible internal branch without computing its old independent bound. Every accepted leaf still checks the requested threshold, and replay checks all splits, spin populations, valence exclusions, and total ionic coverage. There is no weaker acceptance path.

| Method | Nodes | Individual Q leaves | Q configurations in shared leaves | Fully conditioned sources: construction / replay |
|---|---:|---:|---:|---:|
| Previous source-indicator polynomial | 783 | 364 | 16 | 380 / 380 |
| Joint Boolean lines, internal probes retained | 463 | 56 | 324 | 380 / 56 |
| Joint Boolean lines, structural construction | 463 | 56 | 324 | 56 / 56 |
| Joint simplices, structural construction | 335 | 0 | 380 | 0 / 0 |

The last row must be read together with its endpoint count: **168 charts cover400 assignments, including380 ionic and20 valence assignments**. Chart sizes are at most4. Zero singleton leaves and zero fully conditioned transition-source records do not mean zero configuration enumeration.

The final simplex replay counts **862272 metric scalar endpoint evaluations**, **77384 amplitude scalar endpoint evaluations**, and **31764 transition-group evaluations across charts**. The line replay counts668648 metric and58800 amplitude endpoint evaluations, plus its56 singleton leaves. The smaller tree does not establish lower total arithmetic cost or a measured speedup.

The initial simplex endpoint prototype is retained in `coherent_joint_simplex_endpoint_prototype` and explicitly labeled superseded. Its counters are not acceptance evidence for the final compiler.

## The remaining dependence is now measured

The H6 source has1818 CAR terms,1140 coherent transition groups,16 tied metric parameters and30 expanded local charge factors. The metric labels have spatial range2, but distant and multi-site electron transitions unite several such neighborhoods.

Ignoring the separately handled global valence and spin-count conditions:

- **954/1140** transition groups have metric-factor unions spanning all six sites.
- **994/1140** have combined event, spectator, and metric scope spanning all six sites.
- **60** groups have nonconstant affine spectator amplitudes;1080 have constant amplitudes.
- **All1140 grouped amplitudes have a fixed sign** on their allowed fixed-spin source occupations:594 are nonnegative and546 nonpositive, certified by exact rational range bounds.

This removes absolute-value sign ambiguity for this Hamiltonian. It does not make the full weighted-row minimum separable. The scope counts describe this factorization; algebraic cancellation or auxiliary variables might reduce them. They are not lower bounds on the treewidth of every representation.

The next concrete compression target is the factored sum

\[
R(n)=V(n)-\sum_g I_g(n)\,\sigma_g A_g(n)\,Q(n+\Delta_g)
             \prod_{\ell:\,\ell\text{ touches }\Delta_g}
             r_\ell^{f_\ell(q+\Delta_g)-f_\ell(q)},
\]

where exact source-event range bounds certify sigma_g*A_g≥0. A successful elimination experiment must preserve the product factors and the particle-count/valence constraints while evaluating fewer physical endpoint assignments. Its receipt must count intermediate factor entries, elimination operations, and maximum retained boundary scope. Replacing occupation integers with vectors, or hiding the endpoint loop inside a chart, does not satisfy that criterion.

## Energy integration and verification

Both joint-line and joint-simplex gap recipes replace the previous complement proof in the original and hopping-perturbed energy certificates. The energy verifier rebuilds the gap proof from the actual Hamiltonian and requested threshold, and requires the complete valence reference. Upper witnesses and response directions are inherited for this substitution and independently rechecked; they are not newly discovered compressed states.

The original interval width remains **7.416980322924775e−12 Ha** with17 response directions. The perturbed width remains **7.227269910284279e−12 Ha** with18 directions. The perturbation is the previously tested hopping of strength1/50 between the first two localized orbitals, not an additional molecular system. Full energy verification still processes the400-configuration sector.

Focused tests compare every eligible line/simplex chart of an interfering8-mode model against independent CAR actions, including holes, valence endpoints, and amplitude cancellation. Determinant actions are forbidden during chart compilation. Integration checks cover false thresholds, omitted tree coverage, strict selectors, construction restrictions, and node-budget exhaustion.

All **327 marginal regression tests pass in275.466 seconds**. Seven new standard-library certificate replays match every saved replay field. The exact1140-group amplitude-sign receipt is also generated with the standard library. Detailed records are in `results/marginal_final_validation.json` and the accompanying `joint_occupation_progress.json`.

## Remaining blockers

1. **Global proof cost.** The finite gap proof still accounts for every physical endpoint. A genuinely shared elimination or positivity certificate is unproved.
2. **Scalable state and response construction.** The energy witness/response stage retains full-sector work. The valence reference itself grows combinatorially; the current energy implementation has a32-state retained-reference cap.
3. **Transfer across system size and Hamiltonian families.** This result concerns one six-H finite-basis model and one already tested perturbation. It does not establish accuracy/cost scaling on larger, qualitatively different molecules.
4. **General marginal representability.** These weighted operator inequalities certify particular Hamiltonians. They do not characterize the entire cone of physical two-particle marginals.
5. **Materials and synthesis.** Finite-temperature free energies, dynamics, phases, kinetics and realizable synthesis require additional validated machinery beyond a ground-state energy certificate.

The local amplitude/metric decoupling error is solved for the admitted simplex charts. The central global compression problem is still explicit and open.


Follow-up: [Charge-sector spin polynomials](marginal_charge_spin.md) removes the per-spin endpoint loop by sharing exact charge-metric kernels and checking projected positivity. Charge-pattern enumeration remains explicit.
