# Molecular certificates from sparse actions and occupation proof trees

The H4 rectangle now has an independently constructed and replayed interval of width **1.0076048598779602e-10 Hartree** without building a complete sector list or matrix. The input consists only of its rational Hamiltonian, mode count, and particle count. The algorithm starts from one determinant and produces its own retained reference and integer upper witness. It does not import the preceding reference or full-sector eigenvector.

This removes a construction dependency, not the many-body scaling problem. The rectangle's construction evaluates 56 distinct source actions and references 68 of 70 determinants through those actions. Replay evaluates 42 actions and references 60 determinants. The square's construction references every determinant and remains inaccurate under the scalar response approximation.

## The complementary energy proof

For a retained determinant subset P and its coordinate complement Q, real Hermiticity permits the row bound

\[
g_Q(s)=H_{ss}-\sum_{t\notin P,\ t\ne s}|H_{ts}|,
\qquad H_{QQ}\succeq\min_{s\notin P}g_Q(s)I.
\]

The verifier proves a common lower threshold gamma using a binary tree of partial occupation assignments. It never needs a list of all Q determinants. A node is accepted only by a recomputed uniform bound, by an exact terminal row check, or by two verified child branches. Forced occupations from the fixed particle count are resolved exactly. A fully specified retained determinant is excluded explicitly. Every accepted tree covers exactly binomial(M,N) states, including the excluded retained states.

After normal ordering, the diagonal of the admitted two-body Hamiltonian is a quadratic occupation polynomial. At a partial assignment it becomes

\[
d_0+\sum_{i\in U}h_i n_i+\sum_{i<j\in U}J_{ij}n_i n_j,
\qquad\sum_{i\in U}n_i=k.
\]

Its lower bound is d0 plus the k smallest h_i and the k(k-1)/2 smallest J_ij, including zero coefficients for missing pairs. Every possible occupied subset contributes exactly those counts. Minimizing the two sums independently can only underestimate the actual diagonal, even when coefficients are negative.

Each off-diagonal CAR word has an exact activation mask and occupation-flip mask. If its activation is impossible in the node, it contributes nothing. Otherwise its absolute coefficient is included conservatively in the radius, unless **every** possible activated source lands in P. That exception is verified by counting fixed-N completions of the combined node/activation mask and comparing with the number of matching inverse images of the finite retained set under the flip mask. The test does not enumerate the unretained completions.

The triangle inequality makes this safe without relying on cancellation between words. Fully specified complementary leaves instead aggregate actual CAR actions before computing their exact row bound. Tests compare the partial bound against every feasible ascending-mode prefix of both 70-state molecular fixtures using a separate explicit-matrix oracle.

An earlier independent-interval diagonal bound pruned no branches: it required all 44, 53, and 38 complementary rows of the rectangle-coupling, rectangle-witness, and square-coupling references. The fixed-particle-count diagonal bound and exact transition exclusions enable the new pruning.

## Retained lower and upper certificates

Applying H only to retained determinants yields both H_PP and sparse columns W=H_QP. Their exact Gram matrix W^T W requires no additional Q actions. Once the tree certifies H_QQ>=gamma I, the sufficient lower check is

\[
\gamma>b,\qquad
H_{PP}-bI-\frac{W^TW}{\gamma-b}\succ0.
\]

This is the existing scalar Schur certificate, now supplied by sparse construction and a complete tree proof. The independently evaluated upper witness is an integer amplitude dictionary on the current retained determinants. It uses exact CAR actions and a rational Rayleigh quotient. All numerical eigensolves act only on retained matrices of dimension at most 32, and their outputs are proposals rather than accepted proofs.

Selection begins at the lowest-index N-occupied determinant. The counterexample-first variant adds a complementary determinant that fails the current gap threshold. The residual-first variant first follows substantial external Ritz residuals, then resolves gap counterexamples. After obtaining a gap, both variants use the scalar effective matrix to propose additional coupled determinants. These are heuristics; their acceptance depends on the exact lower and upper checks.

## Results and costs

| Fixture and selection | Retained size | Interval width (Ha) | Construction source actions / referenced determinants | Replay source actions / referenced determinants |
|---|---:|---:|---:|---:|
| Rectangle, counterexample first | 28 | 1.0076048598779602e-10 | 56 / 68 | 42 / 60 |
| Rectangle, residual first | 28 | 1.0076048598779602e-10 | 56 / 68 | 42 / 60 |
| Square, residual first | 32 | 0.16418909584125116 | 65 / 70 | 52 / 70 |
| Square, counterexample first | 32 | No accepted certificate | 59 / 70 | — |

A source action count measures distinct determinants on which H is evaluated. The referenced count additionally includes every destination in those sparse action maps. It is the more conservative indicator of configuration-space growth. Avoiding an explicit sector generator is not evidence that the frontier stays small.

The rectangle proof uses 73 tree nodes. Sixteen pruned branches cover 28 Q states; 14 other Q states receive exact row checks. Along with P28 this covers all 70 states. The square's accepted tree uses 103 nodes, 14 pruned branches covering 18 Q states, and 20 exact Q rows.

Initial exports took about 0.36–0.75 seconds for the accepted cases. Recorded independent standard-library replays take about 0.02–0.05 seconds. Construction was repeated from files containing only the Hamiltonian and sector labels; all accepted certificate objects matched exactly. The broader source/destination counts above were recorded in that controlled reconstruction. These are finite-run measurements.

## Validation and boundaries

The three accepted exports have independent `python -S` replays. The square counterexample-first exhaustion is stored separately as `not_accepted`; it is not an energy certificate. The wide square interval is a valid certificate that misses the requested accuracy target.

Tests check sparse actions and exact leakage against full matrices, every feasible prefix bound on both fixtures, complete Q coverage, sparse Rayleigh evaluation, missing branches, false pruning, invalid thresholds and witnesses, changed Hamiltonians, and construction from a bare Hamiltonian with full-matrix calls forbidden. Retained eigensolve sizes are checked during that construction.

A separate 64-mode, 32-particle diagonal sanity test covers its binomial(64,32)-state sector with 65 tree nodes and zero individual Q actions. It verifies that tree coverage is genuinely compressed on a simple family; it does not establish compressed coverage for interacting chemistry.

The implementation has explicit refusal budgets: 32 retained determinants, 4,096 cached source actions, and 100,000 proof nodes, with 2–64 modes admitted. It validates real Hermiticity, number conservation, and a quadratic occupation diagonal. No polynomial bound guarantees that a useful certificate fits those budgets.

The subsequent [optimized residual response](marginal_sparse_response.md) now closes the square precision gap at P32 with12 response directions, width1.00251e-10 Ha, without a full Q matrix or inverse. Its final replay nevertheless evaluates all70 source states. The next construction problem is to retain this directional precision while controlling the action frontier as size grows. The preceding small-block response closes square precision after explicit construction; this tree construction avoids that dependency but currently uses a loose scalar response. A method that simultaneously keeps the construction frontier and response compact remains unproved. General molecular scaling, the asymmetric quartic optimum, and a universal physical-marginal boundary are still open.
