# Eighteen three-spectator constraints integrated into exact energy certificates

All eighteen previously violated three-spectator hopping conditions now enter
the production energy and family verifiers. Both targets have stronger exact
lower certificates. W=0's new lower also strictly exceeds a freshly sharpened
ceiling for the entire preceding two-spectator family. That whole-family
separation has not been established for W=1.

The selected artifacts are in
`results/marginal_graded_hubbard8/three_spectator/<case>/final/`.
`three_spectator/summary.json` records exact fractions, accepted receipts and
remaining limitations; `three_spectator/provenance.json` records source hashes.
The decimals below are rounded displays of exact rational results.

| U=4, t=1, V=1/2 result per site | W=0 | W=1 |
|---|---:|---:|
| New periodic lower | -0.6429317747672409 | -0.6606320334732925 |
| Improvement over preceding accepted lower | 0.000012123576518 | 0.00007916114251224 |
| Enlarged-family ceiling | -0.6429046459320169 | -0.6600346152912886 |
| Certified enlarged-family gap | 0.000027128835223944 | 0.000597418182003887 |
| Positive physical sources in ceiling | 137 | 137 |
| Million-site open lower | -0.6429342747672409 | -0.6606365334732925 |
| Physical variational upper | -0.6106763470511881 | -0.6184244823693281 |

Family ceilings limit the strongest lower certificate in the specified family;
they are not physical ground-energy upper bounds. The physical energy intervals
remain substantially wider. Neither enlarged-family optimum is exactly attained
or tightly resolved beyond the finite gaps displayed here.

## Production extension and preserved gates

For a hopping pair i,j on five sites and the other three positions k,l,m, define
`C = q_k^p q_l^q q_m^r B(i,j)`, with powers in {1,2}. Particle-hole invariance
requires `p+q+r+j-i` odd. Forty primitive products contain four self-reflected
products that cancel; the remaining thirty-six form eighteen reflection pairs.
With `Y=C-reflected(C)`, the six-site correction is `T=Y_left-Y_right`.
Its translated sum vanishes on a periodic chain. The preceding independent
CAR probe established all eighteen directions' independence modulo the old
family, full-Fock symmetry, and a norm bound of eight.

`experiments/marginal_three_spectator_hopping.py` supplies the canonical actions
and bounded exact coefficient parser. ENERGY v16 adds their exact projected
matrix before local PSD testing. FAMILY v12 requires all eighteen moments to
vanish exactly and requires the complete preceding two-spectator hierarchy.
Fixed coefficient fields are forbidden in family inputs, which constrain
moments instead. Older energy versions reject the new field, and the matching
driver rejects preceding-family caps attached to v16 energies.

The new nine-shape family has a source cap of 137; the preceding mode keeps its
119-source cap. Local energy verification still covers all 4096 Fock states in
94 blocks, with maximum local PSD dimension 200. The 4096-character rational
weight limit is unchanged. No physical window, target interaction or opening
cost was added. New correction terms telescope to zero in the periodic proof;
the open-chain conversion still uses the original physical interactions.

Independent reconstruction confirms all eighteen new moments vanish in both
selected mixtures. The thirty two-spectator, fourteen one-spectator, four
pair-transfer, four spin and coherent hopping moments also vanish. Averaged
charge laws retain exact classical Markov extensions. These necessary checks
do not establish an extension of the full quantum density matrices.

## Numerical construction and exact acceptance

Thermal discovery now optimizes 136 coefficients. Each run is bounded by 500
full spectral evaluations at two temperatures, with gradient/curvature checks
and a fresh physical-matrix reconstruction. Initial runs for both targets used
500 evaluations. A subsequent W=0 polish used 390 evaluations and W=1 used 500.
The selected energy recipes are these polished outputs, accepted only after
fresh standard-library PSD replay. Numerical stopping statuses do not establish
optimality. Frozen transfer below uses the initial accepted W=0 recipe.

The first forty-round candidate-feasible trust searches did not produce useful
family caps. Their export then failed because of a missing import in the new
discovery helper. The bug was corrected, failed logs/source versions preserved,
and two global runs affected by the same bug were interrupted and restarted.
This was a construction failure, not an accepted mathematical result. An early
energy replay using obsolete congruence witnesses was likewise interrupted;
all accepted energy replays use freshly generated matching witnesses.

Forty rounds of global pricing, with up to two eigenvectors per physical block,
produced 10514 W=0 candidates and 10794 W=1 candidates. Their numerical LP bases
were reconstructed rationally. W=0's 136-source proposal passed replay with
ceiling -0.6428468597850456. Its longest weight had 4091 characters. W=1's
proposal had a 4217-character weight and the production verifier correctly
refused it. The numerical ceiling -0.6591750977671775 is therefore not an accepted
family certificate.

Further forty-round candidate-feasible pricing passes retained those candidate
ledgers and rounded newly added vectors at integer scale 10^6. Existing vectors
remained intact. The output ledgers have 10869 and 11331 candidates. Exact
reconstruction produced the selected 137-source mixtures. Their longest weight
strings are 2834 and 2964 characters, within the unchanged gate. Both final
production family replays accepted.

Rational reconstruction had become expensive: the global construction runs took
about 387 and 417 seconds. A new discovery-only fraction-free solver clears
column denominators, checks every integer division, then verifies the original
rational equations. It is bounded by 137 rows, 4096-bit cleared inputs and
50000-bit intermediate pivots. It agrees with independent Fraction elimination
on tests, including an overdetermined rational system. The older 85-row helper
and all production acceptance gates remain unchanged. The later construction
runs took about 226 and 207 seconds, but their candidates and pricing strategies
also differ, so this is not a controlled speedup or scalability benchmark.

The last reduced eigenvalues remain negative. Neither global pricing nor trust
pricing has been shown to converge. Exact feasible mixtures certify the stated
ceilings regardless of those numerical search limitations.

## Comparison with the preceding whole family

For W=0, the preceding two-spectator candidate ledger was pruned while preserving
every source in its accepted mixture and all determinant anchors. A bounded
trust search gave a freshly replayed 111-source ceiling of
**-0.6429325359567474**. The new three-spectator lower exceeds this ceiling by
exactly the positive rational recorded in `previous_family_comparison.json`,
approximately **0.00000076118950649 per site**. The comparison checks identical
target, chain size, projector sources, ratio, overlap ceilings and sparse shape
span. This is strict separation from every certificate in that preceding family.

For W=1, the new lower is still below the preceding accepted ceiling
-0.6605062601234778. The signed comparison is -0.000125773349814706. A further
bounded coherent-subspace attempt for the old family returned an inaccurate
numerical objective near -0.66053622; it was not converted into an accepted new
cap and does not establish separation. The new W=1 lower improves on the previous
accepted lower, but improvement over the entire preceding family remains unproved.

## Frozen transfer and ablation

The initial W=0 recipe was transferred to U=5, t=1, V=1/4, W=-1/5, on the same
half-filled chain. Projectors, overlap ceilings, penalties and correction
coefficients were frozen. Physical profiles were rescaled or shifted to match
the target; only the scalar local spectral threshold was recomputed. Three
independent energy replays and a frozen-field comparison were accepted.

| Million-site open-chain recipe | Lower energy per site |
|---|---:|
| New frozen three-spectator recipe | -0.5225841142675748 |
| Same recipe with only eighteen new terms removed | -0.5254907142675748 |
| Previous frozen two-spectator recipe | -0.5226432283437589 |

The physical upper is -0.4885616802989547. The new frozen recipe improves on the
previous frozen recipe by approximately 0.0000591140761841 per site. Removing
only its eighteen new terms loses exactly `14533/5000000 = 0.0029066` per site.
The larger ablation effect reflects a jointly adapted fixed recipe; it is not
a ceiling for a reoptimized family without those terms. This demonstrates finite
coupling transfer, with geometry, filling and interaction range unchanged.

## Validation and outstanding work

All **58 focused production/boundary tests** passed. The full suite passed
**982 tests and 102 subtests in 1099.08 seconds**, with the existing calibration
return-value warning. Three subsequently added fraction-free discovery tests
passed separately in 0.40 seconds. No production source changed after the full
suite was collected; this is not a new full 985-test run.

Selected energy, family, signed comparisons, independent closures and frozen
transfer receipts were audited against current source hashes. Earlier receipts
remain historical evidence, with changed sources checked against preserved
snapshots. Failed, interrupted and nonaccepting construction attempts are kept.
All launched computation and validation jobs are terminal. No GPU, paid resources
or additional agents were used.

The enlarged family gaps remain open, especially for W=1. The new mixtures have
not been tested against every quantum overlap or positivity condition. A useful
next diagnostic is full five-site reduced-density-matrix overlap agreement,
rather than assuming that the tested scalar moments exhaust quantum consistency.
General representability, generic molecular/long-range/higher-dimensional transfer
and computational cost at requested accuracy remain unproved. The goal is active.
