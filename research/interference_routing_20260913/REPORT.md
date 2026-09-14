# Interference-routing investigation — 13 September 2026

**The useful refinement is to certify interference and diagonal slack together, while accounting for every shared positive edge.** I implemented that refinement, reproduced rigorous intervals on three actual H4 geometries, and built a compact routing compiler for a restricted family of sign-frustrated fermion chains. This is a completed bounded investigation. It does not establish a new general classical solution for molecular many-body prediction.

The starting proposal's matrix implication is sound. Its separate requirements that the signed Laplacian be positive and every local energy exceed the endpoint are sufficient but unnecessarily restrictive. Adding a ground vertex turns each diagonal correction into a signed edge. One joint routing matrix can then certify their combined effect, including mixed-sign corrections. The full [derivation and exact counterexamples](/Users/aidenlippert/Documents/Spectra/research/interference_routing_20260913/PROOF.md) distinguish this algebraic extension from any claim about compression.

Two exact four-vertex counterexamples identify unsafe compression shortcuts. First, two negative edges can each pass an isolated resistance test while their shared positive network has a Rayleigh quotient of -3/10. Second, a positive **one-particle** exchange graph can have a **two-particle** eigenvalue of -2. The latter [explicit integer witness](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/sector_counterexample.json) rules out simply replacing a configuration-space proof by a mode-graph proof for occupation-exchange operators. Both examples are minimal in their stated simple-graph/exchange families. Neither contradicts the original proposal's correctly sector-restricted joint matrix inequality; no novelty claim about these examples is made.

Signed-Laplacian/effective-resistance mathematics is established; it is an antecedent, not an invention of this investigation. [Zelazo–Bürger](https://arxiv.org/abs/1408.2187) and [Chen et al.](https://eeqiu.people.ust.hk/wp-content/uploads/2021/09/Characterizing-the-Positive-Semidefiniteness-of-Signed-Laplacians-via-Effective-Resistances.pdf) provide the relevant foundations. The [source audit](/Users/aidenlippert/Documents/Spectra/research/interference_routing_20260913/SOURCES.md) also covers SOS and diagonal compensation.

**The actual failing fixture was the H4 creator-channel certificate.** The current repository already records successful full-dictionary H4–H10 intervals; an older H10 failure status is superseded. I therefore used the still-failing restricted H4 dictionary as the test. Its existing exact dual obstruction excludes a 0.0016-Hartree certificate in that particular coefficient-L1-penalized cone. That obstruction does not exclude a different proof representation.

The unchanged CAR checker was rerun on the restricted certificate. Its Hamiltonian was compared exactly against the frozen straight-chain fixture. Using the same newly evaluated rational upper for the comparison, its interval width is **0.003563810575584488 Ha**. The new routing paths produce:

| Frozen stored Hamiltonian | Routing certificate | Certified width (Ha) | Stored proof bytes |
|---|---|---:|---:|
| Linear H4, 1.4 Å, STO-3G | Guided amplitudes; separated Laplacian test | 4.67282628e-8 | 2,336 |
| Same linear H4 | Uniform lower amplitudes; joint grounded test | 1.40358007e-8 | 2,045 |
| H4 rectangle, 1 by 1.5 Å, STO-3G | Joint grounded test | 1.16385807e-8 | 4,181 |
| H4 square, side 2 Å, STO-3G | Joint grounded test | 1.48492512e-8 | 1,985 |

All four widths are below 0.0016 Ha. The grounded straight-chain bracket is

\[
[-3.6669999700000000,\;-3.6669999559641993\ldots]\;\mathrm{Ha},
\]

with the exact upper fraction in the [receipt](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/molecular/summary.json). Displayed upper decimals are approximate; the checker uses the stored exact fractions.

Uniform lower amplitudes fail the separated signed-Laplacian test on every nontrivial component of linear H4. Joint grounding passes all of them. The square case also accepts two negative diagonal-slack entries, after including their full effect in the joint PSD check. This is a concrete demonstration of the refinement, rather than a claim that the original sufficient theorem was false.

**These molecular results are enumerative oracles.** Every H4 run constructs all 70 determinants. Linear H4 splits into components of sizes 20, 8, 16, 8, 8, 8, 1, 1. Dense component eigensolves propose rational endpoints and the upper vector; exact rational replay alone accepts them. The largest joint routing matrix has dimension 20. Proof bytes exclude the input Hamiltonian and the matrix work reconstructed by the checker. The preceding full-space H4 oracle in the repository already reaches comparable precision. This investigation supplies a different proof representation, not a new best molecular accuracy result or demonstrated speedup over a matched exact solver.

Measured per-case construction, numerical proposal, diagnostics, and exact replay took approximately 0.11–0.23 seconds locally. Exact accepting replay accounted for 0.059–0.094 seconds. These are single-run observations, not a controlled performance comparison. Interpreter startup, artifact writing, prior fixture/integral generation, development, and an earlier uninstrumented H4 inspection are excluded. No expensive original SOS construction is needed to generate the routing proof; the original certificate is used only as the failing comparison. All dense eigensolve dimensions and path work are recorded.

The bounded adapter refuses H6 before constructing its 924-dimensional fixed-N matrix because its explicit enumeration cap is 256 states. That refusal is a resource decision, not evidence that H6 lacks a routing certificate. In general D=binomial(M,N) remains exponential, with dense storage O(D²), component eigensolves/rational elimination cubic in component dimension, and additional rational bit growth. Replacing a full matrix by a tree description does not remove those costs.

**A compact local rule does work when its structural assumptions hold.** For a mode-occupation swap Laplacian K, the exact three-mode inequality is

\[
aK_{12}+bK_{23}-tK_{13}\succeq0
\quad\text{when}\quad a,b,t\ge0,\quad t\le\frac{ab}{a+b}.
\]

The [proof](/Users/aidenlippert/Documents/Spectra/research/interference_routing_20260913/PROOF.md) gives the CAR parity string, the zero cases, and the one- and two-particle blocks. The [compiler](/Users/aidenlippert/Documents/Spectra/research/interference_routing_20260913/local_templates.py) allocates adjacent positive couplings to overlapping triangles. A dynamic program searches five rational route ratios in O(M) work for the fixed grammar. It then evaluates a diagonal-potential lower bound and a uniform fixed-N trial upper with the existing exact O(MN) chain DP. Neither discovery nor accepting replay enumerates the many-body sector.

This family has genuine positive off-diagonal couplings and frustrated sign cycles. However, its uniform control is the known frustrated ferromagnetic J1-J2 chain expressed in occupation pseudospins. It is not newly discovered solvable physics. [Richter et al.](https://arxiv.org/abs/0811.3549) discuss the corresponding one-quarter threshold. The witness language and finite-grammar search are supplied by this implementation; no autonomous learning of that language has been demonstrated.

Twenty-five runs cover five scenarios at 4, 8, 16, 32, and 64 modes. Widths below are in **abstract model energy units, not Hartree**:

| Scenario | 4 modes | 8 | 16 | 32 | 64 |
|---|---:|---:|---:|---:|---:|
| Uniform frustrated control | 0 | 0 | 0 | 0 | 0 |
| Alternating positive capacities | 0 | 0 | 0 | 0 | 0 |
| Added density interaction 0.0001 | 0.00005 | 0.00015 | 0.00035 | 0.00075 | 0.00155 |
| Added density interaction 0.2 | 0.1 | 0.3 | 0.7 | 1.5 | 3.1 |
| Varying one-mode fields | 0.15 | 0.25 | 0.55 | 1.05 | 2.15 |

At 64 modes and 32 particles, the sector contains 1,832,624,140,942,590,534 states. The uniform certificate has 62 routing templates, only 2-by-2 local capacity checks, and a 1,544-byte model-plus-witness file. Its diagonal DP uses 246 states at its widest layer and 15,374 transitions. The route discovery examines 558 transitions, or 734 for the alternating case. Requiring all route ratios to equal one fails the alternating case; the finite search finds a valid allocation. These are structural controls with explicit exact constructions, not molecular performance evidence.

The complete local campaign took **1.432 seconds**, including all 25 discoveries, accepting replays, equal-ratio ablations, and proof-file serialization/writes. Discovery accounted for 0.031 seconds, accepting replay 0.771 seconds, and ablation 0.623 seconds. Interpreter startup and the final summary write are excluded. This accounting rerun produced byte-identical proof files to the preserved initial campaign. Exact records and costs are in the [local campaign](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/local/summary.json). The extensive density-perturbation width has an analytic explanation: delta N(N-1)/M, equal to delta(M/4-1/2) at half filling. Solving sign interference leaves a separate state/energy-variation obstacle.

**Transfer to the supplied molecular representations fails a precise scope test.** Every local template and allowed diagonal term acts on at most three consecutive modes. Exact CAR inspection finds 163, 885, 2,867, and 7,093 canonical terms outside that support for H4, H6, H8, and H10 respectively. The [scope receipts](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/molecular_scope.json) retain explicit nonzero witness terms and fixture hashes. This excludes the implemented template span in the supplied orbital ordering. It does not exclude different orbitals, nonlocal templates, additional residual certificates, or a larger operator-program language.

**Verification is reproducible.** All 29 saved intervals and the sector counterexample passed a fresh standard-library-only process in 1.07 seconds; see [final replay](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/final_replay.json). This uses the same accepting algorithms; independent implementation evidence comes from exact CAR assembly across every particle sector of a nonuniform six-mode model and the four-mode counterexample, exact PSD checks, and independent tree-matrix reconstruction. Fourteen focused tests cover those identities, shared-capacity and sector-lifting counterexamples, zero amplitudes, singular pivots, corrupted endpoints, omitted components, wrong fixtures/sectors, finite-grammar completeness, and extensive error. The five reused chain tests also pass in the NumPy-enabled runtime: [19 passing tests](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/final_tests.log). An initial combined `python -S` invocation failed the existing numerical-spectrum test because NumPy was intentionally unavailable; that log is preserved and is not counted as a pass. The entire repository suite was not run. Existing production solvers and verifiers were not edited.

From `/Users/aidenlippert/Documents/Spectra`:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m unittest research.interference_routing_20260913.test_routing -v
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.interference_routing_20260913.replay_campaign --directory results/interference_routing_20260913 --out /tmp/spectra-routing-replay.json
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.interference_routing_20260913.scope_audit --out /tmp/spectra-routing-scope.json
```

The next consequential mathematical target is a route template that handles an actual nonlocal molecular residual while retaining the full shared-capacity and diagonal-slack proof. A bounded search should start from a witness term in the scope receipt and measure whether it improves the existing molecular interval without determinant enumeration. The current result identifies that missing bridge; it does not assume it exists or quietly replace it with the successful engineered chain.
