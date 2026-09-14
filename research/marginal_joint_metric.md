# Joint metric discovery and projector-free certificates

The latest [bounded-degree coefficient quotient](marginal_bounded_quotient.md) rebuilds the accurate H6 proof and discovers fresh perturbed-H positivity coefficients without physical-configuration evaluation or full-population lifting. An [exact monotone-transfer theorem](marginal_monotone_transfer.md) preserves the −6.264-Ha complement bound throughout the specified hopping interval0≤δ≤0.28551855746986, including mixed-sign targets refused by direct compilation. **360 full tests and nine new independent replays pass.** The inherited metric/direction discovery,400-coordinate H6 quotient, and explicit energy-response work remain uncompressed.

The subsequent [analytic Hubbard construction](marginal_hubbard_polynomial.md) removes configuration enumeration and full-population ideal lifting for a structured chain family. Degree-four positive atoms and degree-two number multipliers give exact zero-residual replay through32 sites. The bound is loose; useful molecular accuracy and compressed energy-response work remain unresolved.

Jointly choosing the metric and its positivity decomposition removes every singleton positive occupation atom from the original and perturbed H6 certificates. A second construction forces the metric to vanish on the valence space and removes the global target projector from symbolic compilation. Both constructions pass exact rational gap and energy replay. These are finite certificates; their discovery still uses all 400 spin-sector configurations.

## What changed

The earlier obstruction in [the polynomial-metric report](marginal_polynomial_metric.md) applies to one fixed metric. Enlarging its positivity cone with all quadratic polynomial squares still fails: an exact signed functional has negative pairing −0.062063581067272874 with the required numerator while its complete 79-dimensional moment matrix is positive semidefinite. Adding all quadratic squares multiplied by D−1 also fails for that fixed metric; a second exact functional has pairing −0.01293456186592471 and both complete moment matrices are positive semidefinite. Each functional also satisfies all 86,509 nonsingleton positive-indicator and 9,909 charge-localizer inequalities, normalization, and the coefficient-residual allowance.

The successful change optimizes the metric and both decompositions together. For a fixed threshold, the compiled numerator depends linearly on the metric coefficients. Consequently this joint search is a linear program in the chosen finite feature and atom spaces. The positivity bound for the metric is fixed to remove its arbitrary scaling. This does not require a bilinear search over independently parametrized metric ratios.

## Exact accepted results

| Construction | Hamiltonian | Complement bound | Metric positivity | Numerator positivity | Singleton positive atoms |
|---|---|---:|---:|---:|---:|
| Joint degree-four metric | Original H6 | −6.24 | 0.000999555127 | 0.019648938817480338 | 0 |
| Same metric, fresh proof | Perturbed H6 | −6.24 | 0.000999555127 | 0.019649205741810692 | 0 |
| Vanishing degree-six metric | Original H6 | −6.264 | 0.000999312001 | 0.01369985018722291 | 0 |
| Same vanishing metric, fresh proof | Perturbed H6 | −6.264 | 0.000999312001 | 0.013701871063934998 | 0 |

All four proofs support the existing energy intervals: width 7.416980322924775e−12 Ha for the original Hamiltonian and 3.9293895949815965e−12 Ha for the perturbation. Each reuses and rechecks the previously constructed upper witness and 17 response directions. These are four new accepted proof routes, not four independently discovered states. Both transfer constructions preserve their source metric and metric-positivity proof exactly and reconstruct a new numerator against the actual perturbed Hamiltonian.

The nonvanishing metric uses 48 reflection-tied degree-four features, rationalized at denominator 10^9. The vanishing metric is constructed as v=D u, with u in the same feature space, reducing charge powers with q_i³=q_i. Its denominator is 2×10^9 and its charge degree is at most six. The vanishing joint search at −6.24 returned a negative numerical optimum; −6.264 supplied a positive proposal and passed exact export. The stronger −6.24 certificates remain separate accepted results.

## Why vanishing on the valence space matters

For q_i=n_iα+n_iβ−1, all valence configurations have q=0. If v(0)=0, then throughout the physical sector

\[
v(q')Q(q')=v(q').
\]

Every transition into the valence space already contributes zero metric weight. The compiler can therefore omit the global target projector exactly. This is tested against all physical rows of a small fixture, including valence sources; target exclusion is retained, not approximated.

For the supported number-conserving Hamiltonians with at most four fermionic operators per term, the projector-free numerator has occupation degree at most d+4 when the metric has occupation degree d. Diagonal terms have degree at most two. Single-excitation events have degree two and grouped amplitudes degree at most one. Double-excitation events have degree four and constant amplitudes. Substitution at changed occupations does not increase metric degree. Fixed-number pruning can lower this upper bound.

On nearest-neighbor Hubbard chains with v=D, the actual numerator degree stays four at 4, 6, 8 and 10 sites. The respective term counts are 73, 260, 507 and 838, with 807, 2,121, 4,067 and 6,645 compiler coefficient products. All four compilations run with determinant actions patched to reject calls. The constant metric instead produces degrees 4, 6 and 8 at 4, 6 and 8 sites. This diagnostic measures compilation, not energy accuracy or positivity at larger sizes. [Exact projector-degree argument](marginal_projector_degree.md).

## What has not compressed

The new H6 proofs still have 1,364 number-multiplier terms in each positivity decomposition. The compiler retains 1,764 numerator monomials. The nonvanishing gap certificate is 785,097 bytes, close to the earlier 790,440-byte proof. Its numerator has 280 positive atoms; 266 cover only three physical configurations each. Avoiding singleton atoms alone is therefore a weak measure of compression.

The vanishing H6 metric requires 10,200,892 compiler coefficient products, compared with 7,992,302 for the nonvanishing metric. Removing the projector eliminates a degree-growth mechanism; it does not make this particular H6 replay faster. Metric richness and expansion cost still matter.

The current exporter completes number identities by lifting monomials to full spin population. This is exact but has sector-sized cost. Discovery uses 400 explicit occupation coordinates and a large indicator dictionary. The full energy replay still accesses 400 spin-sector configurations, 20 retained valence configurations and an explicit response basis. Neither the joint metric nor v(0)=0 removes these costs.

The next targets are a directly generated bounded-degree positivity proof, a controlled response representation without explicit retained configurations, and useful accuracy across larger molecular families. General marginal representability and finite-temperature, dynamical and synthesis predictions remain additional unsolved objectives.

## Evidence

The machine-readable records live under `results/marginal_h6/polynomial_metric/`:

- `joint_metric_proof/` and `joint_metric_transfer/proof/`: exact gap and energy certificates and independent replays.
- `vanishing_metric_proof/` and `vanishing_metric_transfer/proof/`: projector-free counterparts.
- `quadratic_sos_obstruction/` and `localized_quadratic_sos_projected_obstruction/`: exact fixed-old-metric separators and independent replays.
- `projector_scaling_diagnostic.json` and `projector_degree_identity.json`: separately scoped diagnostics.
- `joint_progress.json`: current validation and remaining blockers.

The numerical cubic-square search is unfinished. Its full SDP proposal violated atom constraints; the positive reported objective is not accepted evidence. Eight subsequent valid restricted cutting-plane LPs remained negative, which does not establish a complete cubic-cone obstruction. Narrower joint-cone solver errors and time limits also remain indeterminate.
