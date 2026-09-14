# Spectra cloud continuation checkpoint — 2026-09-12

The user explicitly requested moving this local research project into ChatGPT Work Cloud. This archive contains the project files, not a transferred live Codex process. The existing desktop task is `01a08e58-e30c-7800-8df6-ab713cbd7787`, titled **Identify quantum matter geometry**. The desktop history reader was returning old progress while the local files were newer. Treat files and proof receipts as authoritative.

## Full objective and working rules

Resolve the quantum-marginal certificate blockers: certify or refute the current relaxation's numerical limit, identify and test violated physical positivity constraints, seek compact stronger certificates, and evaluate transfer beyond the matched model. Track remaining unproved scalability and general representability claims explicitly; do not equate finite-model successes with solving general quantum chemistry.

The user's original ambition is an algorithmic chemistry/materials compiler through the geometry of physical reduced density matrices. Continue mathematical and computational attacks, using multiple agents when available. Preserve the full objective. Do not call a finite Hubbard result a universal chemistry solution. Do not restart finished explorations or promote numerical proposals into proofs. Report concrete progress in plain language and keep failed proposals distinguishable from accepted results.

User engineering instructions: preserve correctness and refusal paths, distinguish evidence from assumptions, and use only skills from mattpocock/skills unless the user changes that preference. Relevant global guidance is copied under `migration/agent-instructions/`. No cloud VM purchase is requested; the destination is OpenAI Work Cloud.

## First actions in cloud

1. Verify `migration/TRANSFER_MANIFEST.json` against the extracted files. Use the archive's Spectra directory as the working directory. The original folder was not a Git repository.
2. Read this checkpoint, `research/marginal_six_site_projector.md`, and the latest exact receipt named below. The root README and some older status files are historical and are not the current research entry point.
3. Check whether the cloud environment can execute Python, install the listed scientific dependencies when necessary, and retain/upload updated artifacts. Confirm actual code execution before saying the research has resumed. Do not just repeat this checkpoint.
4. Run the five charged-projector tests, then the full marginal suite. The last completed full suite predates the new charged module. Continue with the coupled overlap attack below.

## Latest fully validated energy milestone

For the half-filled open one-dimensional nearest-neighbor Hubbard model

H = -t sum hopping + U sum double-occupancy + V sum (n_i-1)(n_(i+1)-1),

at U=4, t=1, V=1/2 and N=1,000,000 sites:

* Exact lower/site: `-63102110090479/97656250000000`, approximately -0.646165607326505.
* Physical upper/site: approximately -0.6106763470511881.
* Interval width/site: approximately 0.03548926027531686; 32.4402% narrower than the prior interval.
* Last full regression: **469 tests passed in 347.596 seconds**.

Authoritative files:

* `research/marginal_six_site_projector.md`
* `results/marginal_graded_hubbard8/six_site_projector/refined_certificate.json`
* `results/marginal_graded_hubbard8/six_site_projector/refined_independent_replay.json`
* `results/marginal_graded_hubbard8/six_site_projector_refined_full_validation.log`
* `results/marginal_final_validation.json`

The exact matched replay took 290.226 seconds locally. A descriptive support-size string was corrected after that run; its executed driver snapshot, original hashes, current hashes, and AST-only-string-change explanation are recorded in the receipt. This does not change the arithmetic.

The proof uses a physical cyclic six-site source with 400 integer amplitudes. All six-site profile and even five-site telescoping corrections are capped at density approximately -0.663023831373854 by this source. Five adjacent projectors onto it satisfy sum P_i <= 2621073531/10^9 I. The translated average ceiling is theta_h=0.5242147062. The ten-site Gram has 1,280 columns and maximum spin block dimension180. Local positivity covers all4,096 Fock states in94 reflection/spin sectors, maximum dimension200.

The refined local six-site profile has onsite [a,b,10-a-b,10-a-b,b,a], hopping [p,q,5-2p-2q,q,p], density [d,e,5/2-2d-2e,e,d], with

```
a=.204105 b=3.039489 p=.478826 q=1.260412 d=.147254 e=.372037
half penalty = .213504
accepted local ell = -3.1188936
```

The bulk U4 target requires local onsite mean10/3 and sum20, NOT sum24. Hopping sums5 and density sums5/2. Reflection bases are unnormalized: subtract ell times their column Gram diagonal; normalize numerical matrices by sqrt(column norm products).

## Work completed since that full regression

New module: `experiments/marginal_charged_projectors.py`.
New tests: `tests/test_marginal_charged_projectors.py`.
Five focused tests pass (latest root run0.296 seconds); the new full suite has not yet been run. It should contain474 tests, subject to discovery confirming that count.

`results/marginal_graded_hubbard8/discovery/charged_projector_numeric.py` generates a physical charged source from the active local sector (Nup,Ndown,reflection)=(2,3,+). This has five particles. Spin exchange and staggered particle-hole transforms generate three orthogonal partners in (3,2),(4,3),(3,4). All have equal norm and definite reflection.

The particle-hole map on twelve ordered modes is state s -> 4095 xor s, with phase (-1)^sum_occupied(m+floor(m/2)). Tests independently apply the CAR annihilation sequence on all4,096 determinants. Odd source embedding phases are constant per environment column: Vgraded=Vunsigned*D. They preserve VV-dagger and only conjugate the Gram, so odd source parity does NOT invalidate these projector bounds.

Artifacts under `results/marginal_graded_hubbard8/charged_projector/`:

* `source.json`: integer five-particle source.
* `numeric_proposal.json`: bounded fixed-profile numerical search.
* `certificate.json`: **numerical proposal only** with tentative kind v3; the energy verifier does not yet implement that kind. Do not pass it off as an accepted energy certificate.
* `overlap_replay.json`: exact charged overlap bound accepted in23.475 seconds before the final common-helper refactor. Current focused tests pass after that refactor; replay the actual source again before publishing a new combined result.

Four charged projectors over four windows have1,024 Gram columns, maximum spin block96. The exactly accepted ceiling is B_c=264879339/125000000, giving theta_c=B_c/4=264879339/500000000=0.529758678. No new combined energy improvement has been certified.

## New mathematical result: separate penalties are dominated

This result was derived and audited immediately before migration, and supersedes the prior proposal to optimize independent charged and half-filled penalties.

The five local source vectors occupy disjoint spin/particle sectors, hence P_h+P_c <= I. Suppose

K+kappa P_h+lambda P_c >= ell I,

with nonnegative penalties and separate translated ceilings theta_h and theta_c. Set delta=min(kappa,lambda), reduce both penalties by delta, and replace ell by ell-delta. The new positive matrix equals the old one plus delta(I-P_h-P_c), so it remains positive. Its density bound improves by

delta*(theta_h+theta_c-1)/5.

Here theta_h+theta_c=1.0539733842>1 exactly. Thus every certificate with both penalties positive is dominated by one with one penalty zero, for ANY mean-correct local profile. The charged-only branch cannot cross the old profile ceiling: the cyclic half-filled witness has zero P_c expectation and its K expectation is fixed by the profile means. The independent charged extension therefore cannot improve the best half-only family once that ceiling has already been crossed. This is a structural obstruction, not a failed numerical-search claim.

The small fixed-profile numerical gain (~1e-4/site) used a nearly zero charged penalty and retuned the half penalty. It should not be attributed to charged extendibility.

## Next attack: retain cross-family overlaps

The correct next family is a JOINT bound for P_h+r P_c, not addition of the two separate ceilings. Already implemented, with focused direct-operator tests:

```
joint_overlap_grams(half_source, charged_source, windows)
joint_projector_bound(half_source, charged_source, windows, ratio, ceiling)
```

Both are in `experiments/marginal_charged_projectors.py`. Windows2..4 are admitted. Raw columns contain unnormalized physical vectors. The result retains each column's norm and source index. For positive diagonal D with entries weight/source_norm,

C D C^T <= B I iff B D^(-1) - C^T C >= 0.

This avoids square roots and correctly handles unequal source norms. The joint four-window family has1,280 columns, predicted maximum spin block132; verify those counts rather than assume them. No actual-source joint numerical ceiling or exact replay has yet been computed.

Suggested bounded experiment: compute joint Gram eigenvalue proposals for r=1/2,1,2, with explicit outward margins. For each, numerically optimize nonnegative alpha,beta in

K+(alpha+beta)P_h+beta*r*P_c >= ell I,

using the existing five-window half bound for alpha and the joint four-window bound for beta. The density lower is (ell-alpha*theta_h-beta*theta_joint)/5, minus(2t+abs(V))/N for the open chain. Separate spin sectors still contain at most one rank-one update, preserving the integer rank-one PSD optimization. Implement the production accepting energy path only after a useful proposal, with exact all-Fock and joint-overlap replay and meaningful refusal tests. Do not raise matrix caps to force a result.

## Reproduction commands

Run from the extracted Spectra root:

```sh
OPENBLAS_NUM_THREADS=1 python -S -m unittest tests.test_marginal_charged_projectors -v
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_marginal*.py'
OPENBLAS_NUM_THREADS=1 python -S results/marginal_graded_hubbard8/discovery/six_site_projector.py --certificate refined_certificate.json
```

The current local Python is3.12.2. Scientific package versions are listed separately in `migration/requirements-cloud.txt`. Exact accepting replays use only the standard library under python -S; numerical search and the broader test suite use scientific packages. PySCF and GPU/cuopt belong to older optional investigations and were not installed in the active local environment. Do not assume a GPU is required.

## Historical scope and pitfalls

The original U4,t1,V0 million-site interval remains approximately [-.611636,-.5672624487090691]. Three specified nearest-neighbor U,t,V targets have matched upper/lower receipts in `density_transfer`. Generic molecular interactions, longer range, higher dimensions, thermal/dynamic observables, synthesis, universal representability, and accuracy-versus-cost guarantees remain unresolved.

Source states and algorithms are finite constructions. The million-site upper contracts a specified filtered-block physical state; it is not an exact solution of a million-site wavefunction. The fixed-size transfer contraction does not prove fixed cost as requested accuracy rises.

Agents previously made errors in reflection Gram normalization, full-environment spin counting, derivative indexing, and purported independent tests. Root corrected these. Trust accepted exact receipts only within their scope; inspect independent test construction rather than relying on an agent's summary.

No research compute process was still running when the migration checkpoint was prepared. Do not restart a process merely because an old log says running. The local copy remains intact. Cloud execution is not verified until the receiving task actually opens the archive and runs checks.
