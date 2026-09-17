# H12 acceptance and collective-channel continuation

The best new independently replayed original-model interval is **2.588188188 mHa**, compared with the preserved **12.590172569 mHa**. The 1.6 mHa target is **not met**. No exact family-limit theorem was obtained.

The complete standard-library replay checks the actual MPS upper, singlet and nonsinglet lower witnesses, exact orbital transformation, and its charged coefficient allowance. The Hamiltonian, physical sector and upper are unchanged. This is a continuation from frozen inputs, not a fresh calculation from integrals.

| Complete replay | Width (mHa) | Target met |
|---|---:|---|
| dense_t2_contracted | 2.588188188 | False |
| dense_t2_long | 2.588188188 | False |
| matched_continuation | 2.821202290 | False |

The two 300-second continuations started from the same checkpoint. The four-direction export has an exactly checked width **0.233014102 mHa** smaller than the unchanged-family export. Their selector, construction and replay costs are additional. This is one conditional warm comparison, not evidence of an end-to-end speedup.

**That comparison does not isolate a benefit from the new correlations.** The four-direction run completed 1,408 iterations and the control 1,201. At every positive common logged iteration, the four-direction predicted width was worse. At iteration 1200, the control predicted 2.818957452 mHa and the four-direction run 2.836807406 mHa. The equal-time outcome includes different throughput; its reproducibility has not been established. These intermediate scores are numerical diagnostics, not separately replayed certificates.

The new exact singlet residual allowance is **0.365389793 mHa**, or 14.12% of the complete width. Keep this exported scalar, positive squares, upper and other allowances fixed. The existing norm certificate bounds R between -eta and +eta, so replacing only its scalar lower can improve the result by at most 2 eta. Even that optimistic substitution leaves **1.857408601 mHa**. Residual-only repair of this fixed identity cannot meet 1.6 mHa. This does not limit reoptimization or the full family.

Interpreted in the current orbital labels, the supplied four-term H8 polynomial is already contained in both current H8 and H12 cones. Exact rational coordinates reproduce it in existing highest-weight blocks, and an exact spin-average identity handles its adjoint. The external source bundle was not supplied: its negative dual evaluation, the physical orbital-basis match, and any required transport were not checked. The H12 searches use their own frame-bound inputs.

| Numerical experiment | Best predicted width (mHa) |
|---|---:|
| h12/block_l1 | 4.239038216 |
| h12/control_warmfiles | 4.239038216 |
| h12/dense_t2 | 2.585443077 |
| h12/fixed_scalar | 10.180270000 |
| h12/ideal_l1 | 4.132384406 |
| h12/matched_continuation | 2.818957452 |
| h12/scs_accepted | 3.950302930 |
| h12/scs_control | 3.589937086 |
| h12_cached/dense_t2_8 | 2.638276318 |
| h12_cached/mixed_t2 | 2.758745595 |

The existing solver already selected exports by the residual-penalized lower. The new fixed-Gram ideal LP changed its internal coefficient fit; it improved the warmed control by about 0.107 mHa. The block-scale LP exhausted its time limit and retained its source. A second numerical algorithm and a full accepted-L1 conic formulation were both exercised. Their time-limited outputs did not beat the best candidate and are not optimality proofs.

The selected-channel prototype constructs paired cubic maps through three quartic contractions and updates the normal operator with a Woodbury solve. Four directions add 16 Gram entries; each direction can contain 792 terms. The numerical map agrees with direct polynomial evaluation to approximately 3e-15. The reference accepting path independently expands the actual integer factors.

The later eight-direction and mixed-spin experiments preserved earlier directions. The mixed diagnostic used a 2,520-by-2,520 moment matrix and selected four further directions, each potentially containing 2,520 terms. A small final rank does not erase those discovery or verification costs. No complete global cubic coefficient map was built for these additions, but the existing 320,543 coefficient rows and all original prepared maps remain dependencies.

New metered subprocess stages total **111.00 minutes**, with peak child RSS **2.017 GB**. This includes failed searches, retries and byte-verified input caching after macOS offloading caused read delays. It excludes inherited preparation/discovery, unmetered early unit-test/import time, and lightweight receipt analyses. Those exclusions prevent a complete cold-cost or speedup claim. The prior campaign ledger is linked in the machine-readable accounting. No external compute was purchased. GitHub operations stayed stopped.

All 28 focused regressions passed. These include exact comparisons between the contraction and original CAR algebra, original refusal paths, arbitrary-size integers, an independent literal-ladder oracle on all 256 local occupation states, and a complete two-site/two-electron singlet-basis check. The numerical selector’s negative eigenvalues are not energy gains or exact counterexamples to the present family. The trace diagnostic does not perform exact affine/nullspace repair or exact PSD verification. A stalled search is not an obstruction.

The experimental exact paired-contraction replay produced identical input witnesses, exact endpoints and singlet residual allowance. Its complete time was **670.23 seconds**, against **780.03 seconds** for the original replay (1.164× ratio). This is one replay of each path on this host; it is not a fresh solve or an end-to-end speedup. All unmatched factors, sector checks, orbital allowances and original validations remain on the reference path. The paired contraction is an exactly equivalent component with demonstrated scope, not a solution to the 320,543-row global bookkeeping problem.

The isolated molecular paired component also produced the identical exact polynomial and factor statistics. Original CAR evaluation took **26.786 seconds**; direct integer contraction took **2.357 seconds**, a **11.36×** component ratio. This local ratio must not be applied to the whole verifier. The full replay comparison also includes timing variation in unchanged steps.

The retried numerical singlet-trace diagnostic completed. Its trace normalization/ideal defect was 2.12e-14, but the candidate still acted on the trace kernel at up to 2.84e-05. Even the positive-face-only mixture would weaken the proposed dual energy by 1.193934 mHa, and it does not repair that kernel defect. The suggested 2.359938 mHa family floor is **unverified** and is not an accepted obstruction. Exact affine/nullspace repair and PSD proof remain missing.

A small exact follow-up identified **24 physical singlet-null directions** in the two actual 216-dimensional dictionaries: S_plus a_p_beta and S_plus a_p_alpha_dagger. Their commutators, dictionary coordinates and independence are checked exactly. This identifies part of the trace kernel without relying on numerical eigenvectors. It does not repair the dual, establish which null constraints are forced by the current ideal map, or prove a family obstruction.

The remaining scientific items are H12 accuracy or an exact family obstruction; useful selection at complete measured cost; reducing the retained global bookkeeping and accepting cost; frozen-rule transfer; independent fresh-environment reproduction; and physical-model validation beyond finite-Hamiltonian solver precision. The existing H8/H10, water, orbital-expansion and coupling milestones remain preserved.

See [the mathematics](MATHEMATICS.md), [declared experiments](PROTOCOL.md), [run commands](RUNBOOK.md), and the machine-readable [accounting](../../results/acceptance_channels_20260915/accounting.json).
