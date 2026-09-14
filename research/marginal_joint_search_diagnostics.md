# Diagnosing fresh joint metric discovery

Fresh eight-site gamma=−14 discovery is now exactly accepted. The successful multistage search starts from the bare Hamiltonian, generates atom directions in conditional-moment coordinates, then refines their coefficients in the original coefficient coordinates. It imports no enumerated metric. Both the refined certificate and a smaller crossover certificate pass independent standard-library replay. The earlier investigations below explain how this numerical obstruction was resolved without weakening exact export.

## Direct pricing audit

For each newly priced atom, the discovery driver directly computes its reduced-cost violation from the actual appended atom column and the native equality dual. Every selected atom in the audited runs has the improving sign. Solver primal infeasibility and metric normalization are also recorded. These checks found no pricing-sign or orbit-normalization defect.

The initial normalized coefficient-space optimum uses only a scalar multiple of D. Native simplex previously retained phase-I objective 3.942957 throughout its completed iterations. A native interior-point run moved to 2.957243 and then 2.681700, with nonconstant metric coefficients appearing in later iterations. It still exhausted its 240-second search budget. Moving off the plateau is numerical progress, not certificate acceptance.

Restricting weight-block atoms to positive indicators of degree at most four and charge indicators of degree at most two preserves known feasibility, because the existing weight proof uses this smaller family. The numerator block retains degrees six/four. This controlled run still timed out, with last completed objective 2.464386. The metric remained essentially D-only. Thus removing the higher-degree weight atoms did not independently resolve discovery.

## Conditional physical moments

Let B_j be the sum of occupation monomials in a spin-exchange basis orbit and let S_a be a representative monomial support. Define

T_(a,j) = E_Q[n_(S_a) B_j] / E_Q[n_(S_a)].

All entries are computed from binomial completion counts for the two fixed-spin populations, subtracting valence completions. No physical configurations are listed in constructing T. The implementation checks sampled indicator and charge-localizer columns against direct binomial formulas.

Composing T with the reflection/particle-hole Reynolds map changes the numerical residual objective. Pivoted QR selects 619 rows and checks numerical rank. This rank calculation is only a proposal gate: every exported certificate must pass the original full rational polynomial verifier. Rank loss or a conditioning gate failure is refused, including a small-sector case covered by tests.

In the first eight-site run with these coordinates, the metric changed immediately rather than staying at D. The phase-I objective fell from 44.849958 to 0.630701 across eight completed solves, before the 240-second budget stopped the next solve. Objectives in different coordinates are not comparable numerically. The reduced objective was still nonzero, and no certificate was accepted.

## Exact check of the unfinished metric

The last completed conditional proposal was rounded to rational charge coefficients and independently checked on all 1,106 neutral nonvalence charge patterns, using balanced alternating single-spin runs. Its minimum metric is 108209073/200000000, strictly positive. Its minimum compiled numerator is −43494091/500000000, attained at representative state 58980.

This identifies an actual physical weighted-row violation for that unfinished candidate. Its failed proof search cannot be dismissed as just an export-rounding issue. The pattern enumeration is a diagnostic after algebraic discovery and is not used as a fresh construction step or claimed to scale.

## Continuation support

The discovery driver now checkpoints its selected atom directions with the Hamiltonian digest, gamma, normalization, coordinate choice and block-degree settings. Resume validates masks, degree, fixed-population completion counts, orbit canonicalization and uniqueness. It rebuilds the native LP using those directions; it does not inherit metric coefficients or accepted proof weights. Changed target signatures are refused.

A four-site fresh run and a resumed run both passed full exact export and independent standard-library replay. A changed gamma checkpoint was refused. The earlier timed-out eight-site runs predate these checkpoints; an extended conditional run therefore starts once from the same bare Hamiltonian and saves its directions for subsequent continuation.

Sources: `experiments/marginal_conditional_coordinates.py`, `experiments/marginal_reynolds_pricing.py`, and `results/marginal_graded_hubbard8/discovery/joint_reynolds_search.py`. The search receipts and exact finite diagnostic are retained under `results/marginal_graded_hubbard8/`.

Validation of the production coordinate changes: **377 tests passed in 883.418 seconds**; the focused coordinate/pricing checks passed all six tests. Timings are individual runs under concurrent local work, not controlled speedup measurements.

## Exact fresh discovery and refinement

The extended conditional search generated 4,148 atom directions over 16 completed solves. Its last numerical phase-I objective was 3.976650715269584e-9, but original rational export rejected the candidate. Increasing rationalization precision alone also failed: the weight residual remained about 0.00637 against its 0.0009 floor. Numerical feasibility in the transformed coordinates was not sufficient.

Resuming the same directions in coefficient coordinates with tighter native interior-point tolerances resolved this error in one LP solve, with no additional atom directions. The phase-I objective was 3.2605597890612483e-20. Rational export used metric precision 10^12 and proof denominator 10^14, within the existing verifier limits. Independent `python -S` replay accepted the certificate at `joint_reynolds_coefficient_refinement/proof/`.

The exact weight lower margin is 89999763427/10^14 (approximately 0.000900), and the numerator lower margin is 4992814321/50000000000000 (approximately 0.00009986). Both already subtract the exact coefficient-L1 residual allowance. Search, export and replay contain no physical configuration enumeration. The separate earlier 1,106-pattern diagnostic is not an input to this lineage.

The exporter now permits bounded optional precision while preserving its defaults. In-place exact polynomial accumulation avoids repeatedly copying the growing represented polynomial. A four-site default export remains byte-content equivalent as parsed certificate and receipt to the previous implementation. Ten focused tests passed in 2.512 seconds after these late exporter changes; the 377-test full run preceded them and is not presented as validation of the later revision.

## Crossover compression

A further coefficient solve with native crossover retained all 4,148 search directions but produced a sparser positive decomposition. Full rational export and independent `python -S` replay both accepted it at `joint_reynolds_crossover/proof/`.

| Exact proof terms | Refinement | Crossover |
|---|---:|---:|
| Nonnegative indicators and charge localizers | 19,208 | 6,678 |
| Number-identity terms | 13,875 | 18,248 |
| Combined terms | 33,083 | 24,926 |

This is a 65.2% reduction in nonnegative terms and a 24.7% reduction in combined terms. The number-identity component grows. A separate conservative term-trimming experiment removed only five terms, including constants promoted into the lower bound, and did not provide meaningful compression.

The crossover proof has weight margin 44998795183/50000000000000 and numerator margin 1239247257/12500000000000, approximately 0.000900 and 0.00009914. Its unchanged compiler records zero determinant actions and zero occupation endpoint evaluations, with a degree-six numerator containing 10,017 monomials.

## What this resolves and what remains

Fresh discovery at gamma=−14 is resolved for this Hamiltonian and bounded-degree family. It is a projected-Q lower-bound certificate, not a full ground-energy interval. A Schur lower endpoint must lie below the complement threshold, so gamma=−14 remains too loose to establish an accurate energy near the chain's ground scale. The next mathematical experiment is to locate the stronger-bound limit of this metric family and determine whether a richer metric or operator certificate is necessary.

Certificate growth, larger less symmetric targets, compact energy responses, and general physical-marginal representability remain open. GPU numerical optimization is being prepared as a separate performance experiment; no GPU speedup or cloud execution is claimed by these CPU receipts.

## Subsequent threshold and representation results

The same conditional-discovery/coefficient-refinement procedure also exports a fresh gamma=−12 polynomial certificate. Separately, an exact finite dual excludes gamma=−10 throughout the quadratic charge-metric family. This is a metric-family limit, rather than the earlier numerical export problem.

A change to finite-range positive charge products and joint dynamic-programming verification now gives a much stronger fresh gamma=−4.13 certificate with 68 log parameters, 155 generated charge-pattern cuts, and a roughly 12 KB factor certificate. Its discovery starts from uniform factors and its exact replay does not materialize a complete charge-pattern list. See `research/marginal_joint_charge_dp.md` for the derivation, independent replay, remaining response-rank obstruction and precise scope.


The subsequent local-product/count-profile DP path reaches an exact complement bound of -3.81. A separately verified singlet Schur integration now proves the fixed H8 ground interval [-4.3,-4.20384]; see [the energy report](marginal_singlet_energy.md). These later results supersede earlier frontier statements, while the recorded failed polynomial searches remain valid historical diagnostics.
