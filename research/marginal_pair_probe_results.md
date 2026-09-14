# Targeted pair-transfer test

This experiment tests Path A from the attached proposal: add the ten operators B_ij=T_i T_j-T_i† T_j†, where T_i=a_i† a_{i+5}, to the ten-mode matched model. It also tests a broader pair-transfer span and a cubic control under identical multiplier and solver settings.

## What the hypothesis does and does not say

The earlier data certify the physical ground energy and establish the shortfall of the achieved lower certificates. They do not prove an optimality ceiling for all cubic certificates. We have not constructed a certified dual witness establishing such a ceiling. No pseudo-marginal with the claimed pair-transfer pathology has been exhibited.

A quartic word has four ladder operators; a number-conserving quartic is a two-body operator, not a four-body operator. Multiplying a density by a hopping operator generally gives a quartic, not a cubic. A cubic factor's square can contain degree-six terms, including lower-degree terms after CAR reduction, so cubic positivity is not categorically blind to quartic correlations. N=3 is also odd; odd particle count alone cannot explain the observed difference between N=3 and N=5.

The prior upper witness was approximately 3.28100669557401, not 3.281028. The model experiments do not establish chemical accuracy for every molecule with the same electron count. They also use an arbitrary model energy unit, not automatically Hartree. Prior successes were at three and four particles; mode count must not be read as electron count. Tests use numerical libraries for independent checks, while exact certificate replay does not.

## Algebra of the requested operator

Set A=T_i T_j for distinct pairs i and j. The pair hopping operators are even, act on disjoint modes, and commute. Since T_i^2=0, A^2=0. Hence B=A-A† is anti-Hermitian and

\[
B^\dagger B=A^\dagger A+AA^\dagger.
\]

Writing L_i=i and R_i=i+5 gives

\[
A^\dagger A=n_{R_i}(1-n_{L_i})n_{R_j}(1-n_{L_j}),
\]
\[
AA^\dagger=n_{L_i}(1-n_{R_i})n_{L_j}(1-n_{R_j}).
\]

Thus an individual requested square is diagonal in occupations. Its forward-forward and backward-backward cross terms vanish; it does not itself measure pair-transfer coherence. A joint Gram matrix over different B_ij can still contain nontrivial cross products, so the numerical test uses their entire PSD Gram block rather than dismissing the proposal from the individual-square identity.

## Controlled dictionaries

All runs keep the prior full mixed-cubic dictionaries and every pure triple. They match real Hermitian CAR coefficients through degree eight, restricted to the exactly conserved hopping charges. Number multipliers include all invariant Hermitian words through degree six. This enlargement is common to all three runs and must not be attributed to quartic operators.

- **Control:** cubic factors only, with the enlarged multiplier space.
- **Antisymmetric:** add the ten requested B_ij as one Gram dictionary. A fixed integer transformation maps their polynomial coefficients into a monomial Gram representation for export.
- **Span:** add one 61-dimensional neutral dictionary: the identity, the twenty within-pair one-body words, and forty products T_i^(±) T_j^(±). All cross terms within that block are available.

The span includes the requested operators and can express additional coherence constraints. Neither dictionary is the full quartic SOS hierarchy. A failure of these targeted probes would not rule out other quartic operators; success would not establish a universal topological selection rule.

## Verification contract

The checker now accepts an explicit `operator_degree=4` field. Without it, the original degree-three factor and degree-four multiplier limits remain in force. With it, factors may have degree four and number multipliers degree six. The Hamiltonian remains two-body. Every degree-eight residual coefficient is retained in the exact norm allowance.

The first numerical proposals use SCS with eps=1e-9 and at most 100,000 iterations. After inaccurate outcomes, all three comparisons were run with Clarabel using absolute gap, relative gap, and feasibility tolerances of 1e-9 and at most 200 iterations. Export rounds factors and multipliers rationally. The lower bound is still b minus the full coefficient-l1 residual allowance; numerical success alone is never a certificate.

Five new tests check the nilpotency and diagonal-square identity against independent fermion action, exact reproduction of the proposed polynomial operators, degree opt-in and retention of degree-eight residuals, cancellation of a genuine degree-six number multiplier, and an end-to-end small span discovery/export. All 42 marginal tests passed after correcting boolean ladder labels emitted by the existing adjoint helper at the probe/export boundary.

## Results

All three Clarabel runs reported `optimal`. Their exact exported bounds are:

| Probe | Numerical proposed b | Certified lower bound | Certified interval width |
|---|---:|---:|---:|
| Cubic control | 3.25717841665638 | 3.25715345510872 | 0.02385324046529 |
| Requested antisymmetric operators | 3.25717841685414 | 3.25715845010077 | 0.02384824547324 |
| Broader pair-transfer span | 3.25717840983045 | 3.25715430018711 | 0.02385239538690 |

All share the exact upper witness approximately 3.28100669557401, already known to be within 4.865e-13 of the ground energy. Proposed b values agree within 7e-9. Small differences between the certified lower bounds reflect numerical/export residuals; none approaches closing the 0.02383 separation between the numerical proposals and the physical energy. In particular, a nested dictionary having a slightly lower numerical proposal is not interpreted as a mathematical worsening of the relaxation.

The common multiplier basis has 541 elements and the search matches 1,531 independent invariant Hermitian coefficients through degree eight. The extra Gram dimensions are 10 for the requested operators and 61 for the broader span. This is a sparse, targeted enlargement; neither test includes all quartic charge sectors.

The SCS control and requested-operator runs reported `optimal_inaccurate`. Their exact certified widths were 0.02468764697222 and 0.02471567120027, respectively. They are valid but weaker bounds and are retained as diagnostic artifacts. An earlier SCS span run reached export with pre-fix boolean ladder labels and was rejected by the checker; it produced no accepted certificate. Its result is not used in the comparison. An initial antisymmetric run with the same pre-fix issue was stopped before export and rerun after correction.

Clarabel solve times were 4.89, 5.43, and 8.18 seconds for these individual runs. The SCS control and antisymmetric solve times were approximately 351 and 357 seconds. These are single-run observations under different algorithms, not a general performance claim.

## Conclusion and remaining question

The requested targeted quartic operators did not close the gap in this experiment. Neither did a larger neutral pair-transfer span. The proposed mechanism is therefore not established by the data. There is still no certified upper bound on the optimum of the cubic or quartic relaxation, so this is not an impossibility theorem for either dictionary.

The next informative direction is an exactly feasible dual moment witness, or quartic operators in additional charge sectors rather than only neutral pair transfers. A successful bound there would identify a useful certificate enlargement; it would still not prove that odd particle count causes a universal hierarchy transition.

Five new certificates replayed with exact arithmetic under `python3 -S`. Fourteen earlier certificates from the coefficient, compression, and sector-reference experiments also replayed successfully with the extended checker. The new artifacts are in `results/marginal_pair_probe/`; `manifest.json` hashes the current verifier, proposer, tests, report, and outputs. Earlier manifest hashes remain historical provenance and were not rewritten after source evolution.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_pair_probe --probe antisymmetric --solver CLARABEL
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_pair_probe --probe span --solver CLARABEL
python3 -S -m experiments.marginal_collective results/marginal_pair_probe/span_clarabel.json
```
