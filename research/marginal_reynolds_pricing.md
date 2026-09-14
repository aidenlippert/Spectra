# Symmetry closes the eight-site fixed-metric certificate gap

The imported quadratic charge metric at gamma = −14 now has a complete rational polynomial certificate on the eight-site, half-filled Hubbard chain (U=4, t=1). Independent standard-library replay accepts it. This closes the previously open positivity-decomposition problem for this candidate in the existing degree-six cone. No stronger atom family or verifier relaxation was needed.

This result certifies a lower bound on the specified Q complement, not the full ground-state energy. The metric was originally found by explicit charge-pattern enumeration. Its positivity decomposition was discovered algebraically without enumerating physical configurations. These are different stages with different provenance.

## The equation being certified

The metric v vanishes on the singly occupied valence space and is strictly positive on Q. The exact compiler constructs

K_gamma(s) = (H_ss − gamma)v(s) − sum_(t != s) |H_st| v(t),

with coherent signed CAR amplitudes assembled before absolute values; transitions into valence contribute zero. Positive v and nonnegative K_gamma imply the weighted comparison bound on Q.

Both v and K_gamma receive decompositions into nonnegative occupation indicators, nonnegative (D−1) occupation localizers, and fixed-spin number identities. A coefficient-L1 allowance covers the entire rationalization residual. Acceptance depends on the remaining exact positive lower bounds, not a floating LP status or zero residual.

## What changed in the search

The original spin-exchange coefficient space had 1,988 rows for the degree-six numerator. The target is also invariant under site reflection and particle–hole occupation complementation. For the latter, n_i maps to 1−n_i; in centered charge coordinates q_j=n_(2j)+n_(2j+1)−1 this is simply q_j maps to −q_j. At half filling, D is invariant modulo the number identities.

A Reynolds average over reflection and particle–hole transformations projects into an invariant space with 619 rows. Exact quotient calculations check target invariance. Numerical QR chooses independent rows; it is only a proposal operation. All selected atom weights are expanded through the group, and the original full coefficient verifier decides acceptance. Thus numerical rank errors cannot create an accepted false proof.

The successful fixed-metric run selected 2,330 canonical atom columns over 18 rounds. Numerical feasibility was reached after 66.51 seconds of search; preparation, rational export and part verification brought the numerator run to 132.82 seconds. These are one-run measurements, not complexity bounds. The independently discovered weight proof took 14.50 seconds in its separate run.

The same fixed candidate, positivity floors and degree cap had not completed within 240-second search budgets using either rebuilt simplex masters or a native incremental basis. Conditional physical-Q moment coordinates also timed out. Those numerical failures are retained and do not assert infeasibility.

## Exact acceptance

| Part | Positive indicators | Charge localizers | Number-identity terms | Post-residual lower bound |
|---|---:|---:|---:|---:|
| Weight | 1,186 | 94 | 1,235 | 99999995599 / 1000000000000 |
| Numerator | 3,588 | 1,192 | 9,872 | 499999976999 / 500000000000 |

Weight and numerator residual L1 allowances are respectively 4401/10^12 and 23001/500000000000. The complete proof has 6,060 nonnegative atom terms and 11,107 number-identity terms. It is substantially smaller than the full candidate dictionary but still a large certificate; no compact scaling theorem follows.

The unchanged exact verifier compiled 28 transition groups with zero determinant actions and zero occupation-endpoint evaluations. Separate `python -S` replay completed in 6.37 seconds. Candidate digest: `a080237dbc26e3d4a124389310e1a16287aa8b04ee38473055ba670d05ccff49`.

The old D-only metric has exact weighted-row limit −18 on this chain. The previous implicit certificate reached −18.01. The new −14 result crosses that metric obstruction by changing the metric, while retaining the degree-six positivity cone.

## Remaining experiments and limits

The next construction experiment generates the metric directly from the bare Hamiltonian in the invariant coefficient space. Its centered charge features have definite particle–hole parity; only even features are retained. Both weight and numerator metric blocks are transformed, and exact invariance is checked for each compiled retained column. No imported metric or inherited atom directions are allowed in this experiment.

A separate transfer experiment reuses the complete proof on Hamiltonians with a local symmetry-breaking perturbation. Exact replay, including its full residual allowance, remains the gate. Failure of unchanged-proof reuse does not imply that the perturbed Hamiltonian lacks the bound or that a new proof cannot be found.

Beyond these finite experiments, the unresolved issues remain certificate and search growth, adaptation to less symmetric molecular Hamiltonians, compression of energy witnesses/responses, and general representability. None is established by this lattice Q-bound.

Artifacts: `results/marginal_graded_hubbard8/fixed_reynolds/`, its reproducible `discovery/reynolds_fixed_search.py`, `experiments/marginal_fixed_metric_pricing.py`, and `experiments/marginal_reynolds_pricing.py`. The original metric-discovery provenance is in `quadratic_metric_finite_validation.json` and `charge_metric_diagnostic.json`.

## Fresh discovery outcomes and exact feasibility diagnosis

The reusable symmetry implementation passed two focused tests covering group-expanded columns, both metric blocks, dual pricing and refusal gates. A four-site fresh joint run passed exact export and independent replay in 0.21 seconds of construction.

Both eight-site fresh joint runs at −14 exhausted their 240-second search budgets. The first selected 4,148 atom columns and retained phase-I objective 1.001 at its last completed solve. Hard normalization forbade residuals in the mean-one row; that run selected 4,404 columns but retained objective about 3.94296. Neither result is a certificate or an infeasibility proof.

An exact diagnostic expresses the accepted imported metric in the same 17 independent even charge features. It checks both metric and compiled numerator identities modulo the number equations. The metric's Q mean is 3449999999/3450000000. After normalization, the accepted lower bounds are approximately 0.1 for the metric and 1 for the numerator, comfortably above the joint LP's 0.0009 and 0.0001 floors.

The verifier's small residual allowance does not undermine the feasibility conclusion. For a degree-six residual r=sum c_M n_M, the polynomial r+||r||_1 is a nonnegative combination of constants, n_M and 1−n_M. The latter has the telescoping indicator expansion sum_j (1−n_j) product_(i<j)n_i. These indicators have degree at most six. On the eight-site fixed-spin sector, a singleton would require at least four specifications in each spin, hence at least eight total; all nonzero indicators used here are admitted. Group averaging preserves the invariant target. Thus the full normalized degree-six joint LP has a feasible point. What remains unresolved is efficient discovery of that point from the bare Hamiltonian.

The exact membership, normalization and residual-absorption argument are recorded in `joint_feasibility_diagnostic.json`. They import the accepted metric and do not count as fresh discovery.

## Transfer and compression outcomes

Unchanged-proof replay accepted onsite perturbations ±1/10000 on the left site and Hermitian left-edge hopping perturbations ±1/10000 in both spin channels. These break reflection; the onsite cases also break the relevant particle–hole target invariance. All four have separate accepted `python -S` replay artifacts. Larger positive onsite perturbations 1/1000 and 1/100 failed this unchanged proof's residual-margin check. They do not establish physical failure or an obstruction to adapting the proof.

An exact Shannon-merging pass applied 222 indicator partition identities and passed full replay. It did not reduce the total of 6,060 nonnegative atom terms, so it is retained as a compression diagnostic rather than claimed as certificate-size progress.

Final validation for these implementation changes: **375 tests passed in 387.274 seconds**. Full gap replays remain separate from independent weight-only checks and failed numerical search receipts.

## Fresh discovery update

The gamma=−14 metric is now discovered from the bare Hamiltonian through conditional-moment pricing and coefficient refinement, with no imported enumerated metric. The smaller crossover certificate passes independent exact replay and contains 6,678 nonnegative terms and 18,248 number-identity terms. See `research/marginal_joint_search_diagnostics.md` for lineage, exact margins, validation chronology and remaining threshold/energy blockers. Earlier imported-metric results above retain their original provenance.
