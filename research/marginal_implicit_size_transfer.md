# Implicit energy certificates beyond ten modes

The complete implicit pipeline now accepts a declared even mode count with half filling. At M=12,N=6, repeated targeting gives cycle and mixed-interaction energy widths of **4.28057e-9** and **4.21403e-9**, with final retained dimensions 34 and 64. These extend the earlier dimension-only size probes to validated energy intervals.

The reference remains the matched Hamiltonian at t=1/5. The perturbations break individual pair charges and left/right exchange. This is size transfer within a declared structured-reference family, not discovery of a general chemical reference.

## The complement bound at arbitrary pair count

For m pairs and N=m, a pair-charge sector has d double pairs, d empty pairs, and s=m-2d single pairs. Its single-pair pseudospins decompose into j=s/2,s/2-1,... . In a spin-j block,

\[
H_j(r,r)=\frac{m(m-2)}4+r^2,\qquad
|H_j(r,r+1)|^2=t^2(j-r)(j+r+1).
\]

The diagonal is independent of d. All allowed spins have the same integer or half-integer parity. Embed a nonnegative lowest-energy vector of the j block into j+1 at the same r values, padding the endpoints with zero. Its diagonal contribution is unchanged; the negative hopping magnitudes on shared edges increase. The variational principle gives E_min(j+1)<=E_min(j). Non-strict monotonicity suffices, including the j=0 case.

The reference P is the unique all-single maximal-spin j=m/2 copy. Every complementary block has j<=m/2-1, and the maximal complementary spin is realized in the d=1 maximal-spin block. Therefore the existing rational Jacobi bracket for d=1 supplies a lower bound for the **entire** complement, including the lower-spin d=0 sectors.

Coverage is checked through

\[
\mu(s,j)=\binom{s}{s/2-j}-\binom{s}{s/2-j-1},
\]

and charge-placement multiplicity m!/(d!d!s!). Summing multiplicity times 2j+1 over d,j gives binom(2m,m). Tests verify this identity and exact positivity at the proposed complement threshold for every charge and spin block from m=2 through m=8. The ten-mode threshold is exactly unchanged.

## Generalized fixtures and checks

For m=M/2, the cycle fixture retains centered matched hopping changes s(i-(m-1)/2), cross-pair edge weights -s(1+i/m), and the same local mixed terms with right indices shifted according to m. The original ten-mode Hamiltonians match their saved rational exports exactly.

The M=12 reference embedding and the first cycle space are independently expanded only in tests. Their Hamiltonian actions, Gram matrices, and projected Hamiltonian agree exactly with the canonical Fock oracle. Seven reference columns have exactly zero residual leakage after first enrichment. Production construction and replay use the implicit algebra.

At strength 1/100, the M=12 first-space widths are 8.401692099154303e-5 for the cycle and 9.009460501841786e-5 for the mixed interaction. The cycle's first targeted extension has dimensions 7→17→25 and width 5.652849110201561e-7.

| Fixture | Retained dimensions | Final interval width | Construction seconds |
|---|---|---|---|
| Cycle | 7→17 | 8.401692099154303e-5 | 0.73 |
| Cycle, first target | 7→17→25 | 5.652849110201561e-7 | 122.42 |
| Cycle, second target | 7→17→25→34 | 4.280566327028721e-9 | 201.60 |
| Mixed | 7→46 | 9.009460501841786e-5 | 48.37 |
| Mixed, first target | 7→46→55 | 5.809831029268489e-7 | 75.46 |
| Mixed, second target | 7→46→55→64 | 4.214032940783229e-9 | 149.71 |

Timings come from different implementation stages and include reconstruction for repeated targeting; they are not a controlled scaling benchmark. Final certificate sizes are 14,874 and 19,501 bytes. The respective rational lower endpoints are 53740257427/10000000000 and 53748933163/10000000000. Production replay reconstructs both endpoints without a configuration list.

## Repeated targeting and cost

A certificate may now carry a chain of up to three integer target recipes. Each recipe is interpreted in the preceding reconstructed basis. Replay rebuilds every H0-closed extension, checks its exact orthogonality, and rechecks the final lower and upper endpoints. The legacy single-target format remains supported; tests require identical replay results for both encodings. Ambiguous simultaneous encodings are rejected.

The 64-dimensional retained-space cap remains an explicit computational limit. A failed budget or an unfinished assembly is not an infeasibility result.

Inner products now also group atoms by total right occupation, in addition to their empty/double pattern. Those groups are orthogonal even when H connects them. The original cycle targeted construction took 122.42 seconds; replay using the additional grouping took 9.46 seconds. This is an observed implementation comparison, not an isolated hardware benchmark.

The next optimization replaces pairwise contractions of redundant atoms with exact kernel coordinates within each orthogonal group. For encountered atoms b_i, a growing Gram-LDL decomposition computes coordinates from their combinatorial overlaps. A positive Schur residual adds one direction; zero residual adds none. Newly added orthogonal directions have zero coordinates for every previously cached atom, so old short coordinate vectors stay valid. This is an exact quotient by physical null relations, rather than rank computed from atom dictionaries.

With these coordinates, exact replays take 4.02 seconds for the first cycle target, 34.78 seconds for the second cycle target, 5.04 seconds for the first mixed space, and 20.66 seconds for the first mixed target. Existing saved rational endpoints and pivots agree exactly where an earlier replay is available. Tests explicitly cover cache growth, atom-order independence, physical null relations, CAR action, explicit-oracle matrix agreement, size metadata, and chain-format ambiguity.

The coordinate cache has an explicit 200,000-atom limit. Sector rank and rational bit growth remain costs even without enumerating configurations. The final mixed certificate reaches the current 64-dimensional retained-space limit.

The remaining size problem is the cost and dimension required by further enrichment, and the shrinking margin c-eta-u as the perturbation norm grows. The [full-enrichment error theorem](marginal_enrichment_error_bound.md) controls a complete enrichment chain inside that margin. A separate [selective-targeting theorem](marginal_targeting_convergence.md) now proves conditional energy convergence given a certified excited-state gap; it includes an internal-residual error term for approximate Ritz vectors. Neither theorem guarantees small retained dimension.

The mixed 64-dimensional replay initially finished arithmetic checks but failed while formatting a Schur pivot beyond Python's default 4,300-digit integer conversion limit. A bounded decimal-chunk formatter now permits receipt components up to 65,536 bits without changing the process-wide limit. Tests cover exact large numerator/denominator output, signs, unchanged global settings, and budget rejection. This is an export limit, not failed positivity.

The final mixed 64-dimensional standard-library replay passes in 114.00 seconds; its largest serialized pivot is 10,157 characters. All six new M12 certificates have exact replay receipts, match construction endpoints, and are bound to the declared rational fixtures. The full marginal suite passes 147 tests in 40.183 seconds.

The conditional targeting theorem is now enforced for actual approximate target chains through exact internal-residual and next-vector minimization-error gates. Both final M12 chains have additional convergence exports and independent replays; their original energy endpoints are unchanged. The new recipe certifies an error recurrence, not a bound on the dimensions needed to keep enriching.
