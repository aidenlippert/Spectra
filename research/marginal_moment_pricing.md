# Discovering positivity atoms from coefficient moments

The new restricted-master search recovers the H6 complement certificate using 2,172 active atom columns, selected from generated geometry-only seeds by dual pricing. It passes the existing exact rational verifier. A second implementation removes the need to construct or store the full atom matrix: it prices candidate indicators from a bounded-degree moment map and constructs only selected columns. It also passes exact H6 export from bare Hamiltonian input.

## Restricted master

At a fixed threshold γ, the generated charge metric and its weighted-row numerator depend linearly on 48 metric coefficients. Each positivity block is represented by a growing set of nonnegative atom weights. The initial set contains only generated indicators/localizers whose required support has at most two modes; it contains no active directions or coefficients from an earlier proof.

Artificial signed coefficient residuals make the first LP feasible. The objective minimizes their total absolute value. If y is the equality dual and A_j a currently absent nonnegative atom column with zero objective cost, its reduced cost is −A_jᵀy. Positive A_jᵀy identifies a direction that can improve the restricted master. The implementation adds up to 32 directions per positivity block, then solves again. Residual variables belong only to numerical discovery: the final metric, atom weights and number identities must pass the full rational exporter, including its rigorous residual allowance.

On the complete cached spin-exchange dictionary, 30 rounds reached a zero numerical phase-I objective. The final LP had 421 rows, 3,062 columns including 842 artificial residual columns, and 149,740 nonzeros. It used 2,172 atom directions out of the 67,524 available block directions. Search through the final LP took 29.51 seconds; search plus exact export took 57.83 seconds, excluding the pre-existing cache preparation. This is a smaller LP working set, not elimination of dictionary storage.

## Moment pricing without an atom matrix

Let N(n_M) be the exact fixed-number quotient coordinates of occupation monomial n_M, and P the spin-exchange averaging projection. A dual functional has moments

μ(M) = yᵀ P N(n_M).

For required support R and occupied subset O⊆R, an indicator is

I_(R,O) = n_O ∏_(i∈R\O)(1−n_i).

Its functional value is the finite difference

ℓ(I_(R,O)) = Σ_(T⊆R\O) (−1)^|T| μ(O∪T).

For every fixed R, all occupied patterns can be priced together by the upper-subset Möbius transform. This costs O(|R| 2^|R|) scalar operations for that support, without constructing its coefficient columns.

The charge localizer is L=D−1=−1+Σ_i n_(2i)n_(2i+1). Define

ν(M) = −μ(M) + Σ_i μ(M ∪ {2i,2i+1}).

Applying the same finite difference to ν prices L I_(R,O). Positive indicators require |R|≤6; charge indicators require |R|≤4, so all needed moments have degree at most six. Monomials that exceed either fixed spin population vanish and receive zero moments. The localizer uses double-occupation products, not a linear occupation sum.

The constructor stores the quotient-to-moment map and candidate support/assignment metadata. Only the atom directions selected by pricing are converted into coefficient columns. Combinatorial completion counts exclude empty events and singleton positive indicators. Canonical spin-exchange labels represent averaged pairs. No physical configurations or determinant actions enter this construction.

The direct-moment H6 run used 1,764 occupation moments, a moment map with 44,811 nonzeros, and 48,342 canonical candidate labels. It reached numerical feasibility in 34 restricted-master rounds with 2,268 selected atom columns. Its final LP had 3,158 columns including artificial residuals and 149,291 nonzeros. Preparation took 97.81 seconds, restricted search through the final LP 46.51 seconds, and the whole run including exact export 175.35 seconds. These timings include a fresh metric-basis compilation and moment-map construction, unlike the cached restricted-master timing above.

Full rational export certifies metric lower bound 0.000895336777 and positive numerator lower bound approximately 0.0000890436876683, at complement threshold −6.264. No inherited metric, selected atom directions or full atom matrix were used. The CLI patches full atom preparation and physical-state generation to fail if called.

## Verification and limits

The focused pricing test compares the signs and maximal violations from moment pricing with every explicit candidate atom column on a four-site model, for a deterministic random coefficient dual. It also checks exclusion of already active columns. A two-site test constructs and exactly exports a proof with full atom preparation, determinant actions, physical state generation and full-population lifting disabled. The cached restricted-master test exercises exact export and source/budget refusals. Existing constructor and spin-exchange tests cover the shared feature-generation refactor.

Pricing identities are algebraic; the optimizer and its dual scores are floating-point proposals. Failure to find an improving numerical column is not an exact infeasibility certificate. Candidate labels and moments still grow with the degree and number of modes. The H6 quotient remains the full fixed-spin function space before spin averaging, and the current row budget excludes the configured eight-site run. The separate energy integration still inherits its 400-configuration witness and response construction.

Sources: `experiments/marginal_joint_column_generation.py`, `experiments/marginal_moment_pricing.py`, and shared feature construction in `experiments/marginal_joint_coefficient_constructor.py`. Artifacts: `results/marginal_h6/polynomial_metric/bare_joint/column_generation/` and `moment_pricing/`.

## A separate degree barrier exposed by the scaling check

An exact compilation diagnostic uses a nearest-neighbor Hubbard chain with U=4, t=1 and the metric D[1+q_0²q_1²+q_(m−2)²q_(m−1)²]. This metric is strictly positive on Q. At six sites the compiled numerator has Boolean degree six. At eight sites it has degree eight, with 501 coefficients above degree six whose total absolute value after scaling is 5,440.

Every term admitted by the current positivity proof format has degree at most six, including degree-five number multipliers times a degree-one number identity. Consequently those higher-degree coefficients cannot cancel against that format and remain in its coefficient-L1 residual. This diagnostic is neither a gap certificate nor an infeasibility proof, and it does not rule out other metrics. It shows that simply raising the coordinate budget is insufficient to remove this particular residual: higher-degree proofs or a lower-degree compiled metric family are a separate requirement. See `moment_pricing/degree_scaling_diagnostic.json`.

## Fresh discovery on the perturbed Hamiltonian

The same moment-pricing constructor independently accepts the 1/50 hopping-perturbed H6 Hamiltonian at γ=−6.264. Its input contains only the target Hamiltonian and sector. This run builds its own metric features, moment map and selected atom columns; it does not reuse the original solution or invoke monotone transfer.

It reaches feasibility in 30 rounds with 2,140 selected atom columns, 3,030 total master columns and 140,389 nonzeros. Preparation takes 113.58 seconds, search through the final LP 32.43 seconds, and preparation/search/export together 176.99 seconds. The exact exported metric lower bound is 0.000888141702 and the numerator lower bound approximately 0.0000604814667831. These are separate finite-input successes, not evidence of molecular-size scaling.

The three newly constructed gap proofs (cached restricted master, direct moments, and fresh perturbed direct moments) independently replay with the standard library. Each also integrates with the existing energy certificate. The original interval width remains 7.416980322924775e−12; the perturbed interval width remains 3.9293895949815965e−12. The energy reference, responses and 400-configuration witness are inherited in these integrations. All 368 marginal tests pass.

A further exact basis-count diagnostic gives a concrete next step. The current tensor quotient has 4,900 coordinates at eight sites. Restricting its free basis monomials to total degree at most six leaves 3,920, below the current 4,096-coordinate budget. The analogous counts are 400→400 at six sites and 63,504→23,340 at ten sites. This restriction is appropriate only when every input polynomial lies within the stated total degree; the degree-eight diagnostic above cannot simply be projected away. No larger-system positivity solve is claimed by these counts.
