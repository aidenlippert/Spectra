# Connected dimer moments: fixed-order compression and its energy limit

For the uniform open half-filled Hubbard chain, the implemented compiler computes exact normalized moments through order eight of one product-of-valence-singlets boundary using clusters of at most ten sites. The same five clusters serve twelve sites or a million sites. Their total is 2,045 source actions in a local Clifford frame, with 1,306 in the largest cluster. No whole-chain Hamiltonian or configuration list is constructed in this route.

This is a fixed-order result for a declared model and boundary. It does not compute a million-site ground state, certify general physical marginals, or establish efficiency as moment order increases. The extension to arbitrary even chain length uses the analytic lemma below; the code checks finite contractions and identities, not a machine-checked proof of that lemma.

The earlier eight-site ground-energy certificate is unchanged. Its strongest interval has lower endpoint -4.235845 and exact upper endpoint approximately -4.235805939245279. It still uses complete singlet information and a separately certified complement bound, described in [the symmetry-moment report](marginal_symmetry_moments.md).

## Exact construction and proof

Write N=2B and use the normalized boundary

\[
|\Omega_B\rangle=\bigotimes_{j=1}^B
\frac{|1001\rangle-|0110\rangle}{\sqrt2}.
\]

Each four-mode block has even parity and exactly two particles. The model is the open nearest-neighbor chain with uniform real hopping t and onsite repulsion U. The current API requires even N, nonnegative exact rational U,t, and order at most eight.

Define m_n(B)=⟨Ω_B|H_B^n|Ω_B⟩ and formal cumulants by

\[
\log\langle\Omega_B|e^{zH_B}|\Omega_B\rangle
=\sum_{n\ge1}\kappa_n(B)\frac{z^n}{n!}.
\]

No convergence or analytic continuation in z is needed for the finite formal coefficients.

**Connected interval lemma.** The order-n cumulant receives contributions only from intervals containing at most floor(n/2)+1 dimers.

1. Give each inter-dimer hopping edge its own formal coupling. Setting an edge to zero splits the Hamiltonian into even, commuting operators on disjoint components. The product boundary factorizes. Its exponential expectation therefore factorizes, and its formal logarithm is additive. A logarithmic coefficient involving multiple edges can survive only on a connected interval.
2. In a matrix element of an ordered hopping word, the initial and final charge on each side of any cut are equal. Every used inter-dimer cut must consequently be crossed equally often in both directions, at least twice in total. Onsite terms and intra-dimer hopping do not change that cut charge. Products of moments in the cumulant recurrence retain this requirement.
3. An interval of b dimers uses b-1 distinct cuts. It therefore needs n≥2(b-1). Noncommutativity within an interval causes no problem: the factorization in step 1 is used only between disjoint even components.

Translation invariance on the fixed dimer partition makes an interval's weight independent of position. With κ_n(0)=0, extract its exact weight by the finite difference

\[
w_n(b)=\kappa_n(b)-2\kappa_n(b-1)+\kappa_n(b-2),
\]

where w_n(1)=κ_n(1). Then

\[
\kappa_n(B)=\sum_{b=1}^{\min(B,\lfloor n/2\rfloor+1)}(B-b+1)w_n(b).
\]

The moment/cumulant recurrences are exact integer or rational arithmetic:

\[
\kappa_n=m_n-\sum_{j=1}^{n-1}{n-1\choose j-1}\kappa_jm_{n-j},\qquad
m_n=\sum_{j=1}^{n}{n-1\choose j-1}\kappa_jm_{n-j},\quad m_0=1.
\]

Order eight therefore requires only b≤5, or ten sites. Arithmetic operation counts after cluster contraction are independent of N at fixed order; integer bit lengths still grow with log N. The present input bound is N≤10^9. General graphs, long-range hopping, non-product boundaries, and other charge-changing interactions require a new connectivity analysis.

## Why the Clifford frame helps

An exact local Clifford circuit maps the frame vacuum to each normalized valence singlet. Jordan–Wigner conversion uses raw real Pauli words X^x Z^z and exact rational coefficients. Hamiltonian action in this frame starts the product boundary at one coordinate, rather than 2^B physical determinant amplitudes.

The frame retains fixed limits of 4,096 source actions and 65,536 vector coordinates. It does not remove general state growth: a whole twelve-site eighth moment passes with 2,557 source actions, while the sixteen-site eighth-moment probe refuses at the source cap. Connected intervals remove that whole-chain expansion for the declared scalar moments.

The exact frame also passes four-site comparisons for a Hamiltonian with spin-flip and pair-transfer terms that break the older symmetry gates. That checks the change of basis only. Neither the interval lemma nor the existing Lieb-qualified energy certificate is thereby transferred to those models.

## What the resulting energy experiment says

Moments through order seven evaluate any degree-three polynomial Rayleigh witness exactly:

\[
E_0\le\frac{\sum_{i,j=0}^{3}c_ic_jm_{i+j+1}}
{\sum_{i,j=0}^{3}c_ic_jm_{i+j}},
\]

provided its norm is positive. Numerical discovery supplies four coefficients; replay uses their exact rational values and recomputes every moment. The resulting upper bounds are:

| Sites | Polynomial upper energy per site, approximately |
|---:|---:|
| 8 | -0.4992308051 |
| 12 | -0.4925771182 |
| 64 | -0.3428056978 |
| 1,000 | -0.1090593240 |
| 1,000,000 | -0.0036833987 |

Every displayed number is rounded from a rational witness, not a claimed ground energy. The million-site moment contraction plus four-dimensional numerical proposal took about 0.2 seconds in recorded CPU runs. These are measurements of this fixed-order task, not equal-accuracy ground-state benchmarks.

**The deterioration has an analytic explanation.** For U=4,t=1, κ_1(B)=0 and κ_2(B)=5B-1. At every fixed order the connected formula makes κ_n(B)=O(B). Thus the moments of H/√B through any fixed order tend to those of a centered Gaussian with variance five. For the degree-three Krylov space, the limiting Gram matrix is positive definite, so its lowest generalized eigenvalue converges to the smallest zero of the fourth probabilists' Hermite polynomial, scaled by √5:

\[
E_{\mathrm{Krylov},3}(B)
=-\sqrt{5(3+\sqrt6)}\sqrt B+o(\sqrt B).
\]

Its energy per site approaches zero. In contrast, a product of exact two-site Hubbard ground states is a legal fixed-charge trial state. Inter-dimer hopping has zero expectation and

\[
E_0(2B)\le B\frac{U-\sqrt{U^2+16t^2}}2
=(2-2\sqrt2)B\quad (U=4,t=1).
\]

That trial energy per site is 1-√2≈-0.41421356, independent of size. Therefore the optimally chosen fixed-degree scalar Krylov space itself fails to retain energy-density accuracy; coefficient optimization cannot cure this failure. A product of local correlated states escapes the obstruction because it is not a fixed-degree polynomial of the total Hamiltonian applied to the original boundary.

The driver also constructs an explicit rational version of this product witness: use its two-site polynomial coefficients on every dimer. Replay recomputes the local quotient and multiplies it by B; each factor has fixed charge two, so cross-dimer hopping expectations vanish exactly. Its upper energy per site is approximately -0.414213562373095 throughout the recorded sizes. This certifies an extensive trial energy with constant local data, while leaving the interaction error and a matching lower bound unresolved.

A finite moment matrix also supplies no global ground-energy lower bound on its own. A boundary can miss an orthogonal lower-energy subspace, and finite moments do not uniquely determine even its visible spectral support. In particular, the smallest Hankel generalized eigenvalue is a variational upper bound, not a certified lower endpoint. The H8 complement and Schur certificates contain additional global operator information that these scalar moments lack.

## Verification and remaining work

The focused tests check exact original-CAR contraction on two, four, and six sites, including rational parameters and an odd number of dimers; twelve-site agreement with direct Clifford contraction; normalization; zero hopping; moment/cumulant inversion; fixed cluster work at large N; and input/budget refusal paths. The four-site Clifford tests exhaust all 256 frame coordinates for both the Hubbard and modified Hamiltonian.

Final regression validation: **404 tests passed in 288.586 seconds**. The ten focused standard-library tests passed in 1.508 seconds. Logs are `connected_dimer_full_validation.log` and `connected_dimer_focused_validation.log` under the H8 result directory. The subsequent addition of the product-witness fields affects only the research driver and is covered by its separate standard-library replay.

The research driver records eleven chains from two to a million sites. Separate standard-library replay recomputes all nine moments, positive norms, exact Rayleigh quotients, and cluster work. Finite checks corroborate the implementation; the analytic lemma supplies the arbitrary-size extension.

The next unresolved requirement is a compact representation of correlations throughout the system together with a global positive-operator certificate. Candidate local or tensor representations must be tested at equal energy-density accuracy as size grows. They must also supply a controlled truncation bound and survive loss of the matched model's exact symmetries. General N-representability and molecular transfer remain open.

Follow-up: [local energy-chain certificates](marginal_local_energy_chain.md) now supply a coarse two-sided global interval using exact all-Fock local positivity and overlapping windows. The million-site width with six-site certificates is about0.1822 per site. This resolves existence of a fixed-cost coarse global certificate for this chain; sharpness and efficient dependence on requested accuracy remain unresolved.

A separate exact higher-order probe finds w_n(4)=0 through n=9 and w_10(4)=26495 at U=4,t=1. This is additional cancellation on the eight-site cluster, not permission to omit larger clusters at higher orders. Its recomputation is retained in `connected_dimer_moments/higher_order_diagnostic.json`; the production order cap remains eight.

Cluster and cumulant methods are established tools; this report claims an implemented, explicitly scoped instance and its checked consequences, not invention of the method. For related rigorous short-time cluster algorithms, see [Wild and Alhambra, PRX Quantum 4, 020340 (2023)](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.020340). Their short-time result does not establish a general ground-state compiler.

## Reproduction

```sh
python results/marginal_graded_hubbard8/discovery/connected_dimer_probe.py
python -S results/marginal_graded_hubbard8/discovery/connected_dimer_probe.py --replay
python -S -m unittest tests.test_marginal_connected_dimer_moments tests.test_marginal_clifford_moments -v
```

Implementations: `experiments/marginal_clifford_moments.py` and `experiments/marginal_connected_dimer_moments.py`. Data: `results/marginal_graded_hubbard8/connected_dimer_moments/receipt.json` and `stdlib_replay.json`. Earlier whole-chain frame probes remain in `clifford_moments/scaling.json` and `clifford_moments/h8_comparison.json`.
