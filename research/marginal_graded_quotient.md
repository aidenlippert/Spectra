# Total-degree quotients and an eight-site certificate

The joint coefficient quotient now respects total polynomial degree. This removes unnecessary tensor coordinates and, with degree-compatible metric features, enables fresh moment-priced certificate discovery on an eight-site Hubbard chain. The existing rational verifier accepts complement bounds −19 and −18.01. These are larger structured-lattice gap certificates, not eight-hydrogen molecular calculations or new energy intervals.

## The graded quotient

For a one-spin fixed-number slice, the quotient filtration through degree k has dimension C(m,k), for k within the supported particle/hole range. Its free basis therefore has

h_k = C(m,k)−C(m,k−1)

basis elements of degree k, with h_0=1. The exact elimination order preserves polynomial degree. Tensoring the two spin bases and retaining only elements of total degree at most d gives

Σ_(i+j≤d) h_i h_j

coordinates. They span all reductions of input polynomials of total degree at most d and remain independent as a subset of the full tensor basis. Thus no information in that input space is discarded.

`CoefficientQuotient` now validates total degree before reduction, checks the predicted rank, and refuses out-of-degree input. A degree-eight numerator cannot be silently projected into a degree-six quotient. This is a conservative interface: even a special high-degree polynomial that could simplify further is refused at this boundary.

At eight sites, the full tensor quotient has 4,900 coordinates. Total degree six needs 3,920; total degree four needs 1,260. Spin exchange further reduces the degree-four coordinates to 644 per positivity block. H6 uses the same 400 coordinates as before at total degree six.

## Matching metric and proof degree

Feature generation now exposes the charge-feature degree separately from the coefficient/positivity degree. The moment constructor generates positive indicators through the chosen proof degree and charge localizers through two degrees less. It constructs the moment map sparsely rather than allocating a dense quotient-by-monomial rectangle. The final proof still uses the unchanged exact polynomial verifier, whose maximum degree remains six.

The first eight-site run used D times quadratic charge features and degree-six proofs at γ=−16. It reached the 180-second search limit without feasibility; its recorded phase-I objective remained 1.001 while directions were still improving the dual. This is a numerical/resource outcome, not an infeasibility proof. That model has 3,920 quotient coordinates, 13,829 moments, 876,325 moment-map nonzeros, and 358,314 candidate labels.

The smaller search uses only a positive scalar multiple of D. Its numerator has degree four for the nearest-neighbor Hubbard Hamiltonian. The model has 1,260 quotient coordinates, 2,517 moments, 39,594 moment-map nonzeros and 17,386 candidate labels. It never materializes the full positivity matrix.

For U=4, t=1 on eight sites:

| Requested complement bound | Selected atom columns | Search rounds | Whole construction including export |
|---|---:|---:|---:|
| −19 | 2,948 | 39 | 19.05 seconds |
| −18.01 | 3,588 | 46 | 27.37 seconds |

The earlier direct analytic bound for this same parameter choice is −24.001. Both new results pass full rational export and independent standard-library replay. Metric generation and atom selection start from the bare Hamiltonian, with physical state generation and full atom preparation disabled. Timing is from individual recorded runs, not a general complexity guarantee.

## What the fixed metric can and cannot achieve

An independent integer dynamic program determines the exact D-weighted row limit. It tracks prefix alpha/beta populations, double occupancy and the last local occupation. For fixed total double occupancy D, each allowed nearest-neighbor hop contributes D plus its change in double occupancy. A hop into the valence space contributes zero because its target metric vanishes. Maximizing these local contributions at fixed populations gives the worst weighted hopping row without enumerating complete configurations.

The exact row minima for D=1,2,3,4 are respectively −18, −5, 5/3 and 11/2. One critical occupation pattern is [double, up, down, up, down, up, down, empty]. Its diagonal energy is 4 and its D-weighted hopping sum is 22. At γ=−17.99 the exact compiled numerator on this physical state is −0.01.

Consequently no nonnegative numerator decomposition can certify −17.99 with this fixed D metric or a positive scalar multiple. This obstruction applies regardless of how many positivity atoms are added. It does not exclude another metric or a stronger operator bound. The actual complement eigenvalue is not determined by this row calculation. The accepted −18.01 certificate is within 0.01 of the fixed-metric weighted-row limit, but this alone does not prove attainment of the endpoint in the restricted degree-four cone.

The diagnostic recurrence is cross-checked against every four-site fixed-spin polynomial row, and its eight-site critical row against the exact compiled polynomial. Its script and rational receipt are retained under `results/marginal_graded_hubbard8/`.

## Remaining work

A useful bound above −18 on this eight-site example requires changing the metric or the operator certificate. The earlier quadratic-feature search at −16 timed out; the subsequent fixed quadratic candidate now has an exact −14 certificate, as described in `research/marginal_reynolds_pricing.md`. Molecular Hamiltonians can require higher numerator degree than this chain, so the degree-four result is not a general solution to their degree barrier. Candidate scanning and coefficient dimensions still grow. No larger energy witness or response compression is supplied here; the existing H6 energy path still uses 400 configurations.

Sources: `experiments/marginal_joint_coefficient_constructor.py`, `experiments/marginal_moment_pricing.py`. New focused coverage in `tests/test_marginal_graded_quotient.py` checks the exact rank, small-sector polynomial equivalence, fixed-number identities, total-degree refusal, and degree-compatible eight-site feature construction. Existing constructor, pricing and symmetry tests also pass their focused run.

## Distinguishing the metric family from the certificate search

Seeding the quadratic-metric search at γ=−16 with 635 canonical atom directions from the new −18.01 certificate did not produce a proof within its 120-second search budget. It solved fresh metric/positivity weights but inherited those seed directions, so this was not an unseeded discovery run. The outcome remains a solver/resource result.

A separate finite charge-pattern calculation shows that the quadratic metric family is physically capable of a stronger bound. For the checked uniform Hubbard chain, all 1,106 neutral nonvalence charge patterns were enumerated. Single-occupation runs can each alternate spins; their odd-run count is even, allowing global Sz=0. Since target metrics are positive, this simultaneously maximizes the hopping penalty on all single-single edges. The other allowed hopping counts depend only on the charge pattern. Thus one balanced alternating representative per pattern realizes a worst spin row.

A numerical metric proposal at γ=−14 was rounded to exact rational charge coefficients. The actual Hamiltonian compiler then evaluated every representative with rational arithmetic. The minimum metric is 21367273/125000000 and the minimum numerator is 1340465967/500000000, both strictly positive. The source Hamiltonian is explicitly checked against the same U=4,t=1 chain. The candidate and reproducible validation script are retained as `quadratic_metric_candidate.json` and `discovery/validate_charge_metric.py`.

This finite weighted-row validation was obtained with explicit pattern enumeration; it was not itself an implicit polynomial certificate. Subsequent reflection/particle-hole moment pricing has now produced a degree-six decomposition for this exact candidate, accepted by independent rational replay. See `research/marginal_reynolds_pricing.md`. The original discovery still has enumeration provenance; the later certificate does not change that history.

Historical validation at the graded-quotient milestone: 371 full marginal tests passed. Its two polynomial gap certificates passed independent standard-library replay. The fixed-metric row limit and the −14 finite-pattern check are recorded separately as exact diagnostics, not as additional implicit-certificate replays.

## Fresh discovery update

The gamma=−14 metric is now discovered from the bare Hamiltonian through conditional-moment pricing and coefficient refinement, with no imported enumerated metric. The smaller crossover certificate passes independent exact replay and contains 6,678 nonnegative terms and 18,248 number-identity terms. See `research/marginal_joint_search_diagnostics.md` for lineage, exact margins, validation chronology and remaining threshold/energy blockers. Earlier imported-metric results above retain their original provenance.
