# V4: learned experimental planning did not compound

The second learned policy update supplied **zero additional risk reduction** over the first retained rule on all 32 held-out contexts. The exact signed-moment planner attained the optimum in every context. This is a useful null result, not proof of compounding scientific capability.

## What ran

A deterministic CART regressor learned two-shot planning-value corrections from 24 solved design cases. A second tree learned residual corrections from 32 further cases. Training systems had two or three qubits; evaluation systems had three or four, with distinct seeds, noise values and gate combinations. A supplied finite candidate set contains eight signed coupled-Pauli mechanisms. Exact rational Born probabilities describe noisy measurements; two simulated calibration outcomes update the prior in every context.

The learned object is an action-scoring function. The physical grammar, likelihoods and calibrated noise are supplied; the two training sets are not two independently discovered physical mechanisms. Teacher labels come from an exact computational planner. This is an exploratory split, not an externally blinded or preregistered test. The 32 contexts establish exact conditional risks for these cases, not a population confidence bound.

Each proposed first action receives an optimal final-shot continuation. A separate checker integrates its four terminal joint outcomes and verifies the claimed risk without invoking the policy generator. The stronger exact baseline computes signed likelihood moments, selects both actions directly and compiles a certificate without redundant continuation search.

## Results

Lower risk is better. Times below are measured implementation costs, not mathematical lower bounds or hardware benchmarks.

| Method | Mean exact error | Mean regret versus optimum | Mean planning/verification time |
|---|---:|---:|---:|
| One-step goal heuristic | 0.402593 | 0.011187 | 2.186 ms |
| Fixed heuristic portfolio | 0.397481 | 0.006074 | 2.985 ms |
| First learned rule, retained | 0.395526 | 0.004119 | 3.618 ms |
| Two cumulative learned rules | 0.395526 | 0.004119 | 3.909 ms |
| Second-stage learning from scratch | 0.397027 | 0.005621 | 3.830 ms |
| Conventional batch learner | 0.395315 | 0.003908 | 3.799 ms |
| Nearest-record retrieval | 0.395474 | 0.004068 | 3.646 ms |
| Exact signed-moment planner | 0.391406 | 0.000000 | 1.941 ms |

The raw second-stage score was slightly worse than the raw first-stage score: mean error rose from 0.398733 to 0.399092. Adding the second score to a verified candidate portfolio preserved the first score's result but did not improve it. A conventional batch learner and nearest-record retrieval were slightly better than cumulative learning. Exact context retrieval had zero cache hits. Feature erasure and shuffled-label controls are included in the [complete summary](summary.json); no control supports a compounding claim.

## Resource ledger

Training consumed 112 simulated calibration shots and 280 qubit preparations across 56 contexts. Teacher calculation used 4,032 weighted-likelihood products and 20,384 signed-Gram products, taking approximately 0.0648 seconds in this run. The first and incremental tree fits took approximately 0.0546 and 0.0655 seconds. Other baselines' training costs are recorded separately.

The 32 evaluation contexts consumed another 64 simulated calibration shots and 224 qubit preparations. The same context record is supplied to every method. A selected policy **would use** two target shots, 2n qubit preparations and the recorded one/two-gate source executions. Target-policy error was integrated exactly; those target shots were not actually sampled. No hardware experiment was performed.

Per-policy records include candidate count, likelihood multiplications, posterior divisions, rational operand sizes, teacher work, prediction comparisons/distance coordinates and independent verification work. Feature extraction is charged to the heuristic/learned paths. Source-table construction and calibration are common setup; their full bit-operation and physical implementation costs are not benchmarked. CART split counters are algorithm events, not complete floating-point operation counts. There is no resource-complete scientific advantage claim and no favorable amortization threshold: the learned cumulative implementation was slower at reuse and less accurate than the exact baseline.

## Exact obstruction discovered during the attack

Consider hidden independent fair bits u,v,w and target sign (-1)^(u+v). Four actions return (u,v,w,w) in one world and (w,w,u,v) in another. All twelve per-action feature vectors used by the learner are identical across these worlds. Yet their exact two-shot first-action risk vectors are respectively

`(0, 0, 1/2, 1/2)` and `(1/2, 1/2, 0, 0)`.

Every rule restricted to this feature interface incurs mean regret 1/4 over the two worlds. This is an exact impossibility result for that interface, including arbitrary training and randomized decisions. Cross-action information repairs the ambiguity. It is not a proof that this larger deterministic-channel construction alone explains the narrower benchmark's null.

The [signed-moment theorem](../../research/v4/signed_moment_theorem.md) gives the exact sufficient information for two-shot binary prediction. The [higher-order theorem](../../research/v4/higher_order_obstruction.md) proves that no fixed interaction order handles all growing horizons; it also states a charged labelled-sample estimation bound. These are elementary finite-model derivations, with [verified prior-art context](../../research/v4/moment_sources.md), not claims of new physical law or unrestricted novelty.

## What remains open

The next positive result must acquire relevant model or representation information from observations, use it to alter a subsequent discovery procedure, and improve at least two successive held-out stages after acquisition and verification costs. A stateful conventional learner is allowed to tie; beating an emulator is not a meaningful requirement. Equality with a frozen first-stage learner, as observed here, does not pass.

The broad goal remains active. V2's positive complementary-acquisition theorem survives; V3's strong-damping approximation remains valid. Neither this failed policy learner nor the moment formula supplies the missing evidence of learned scientific structure in general coupled matter.

Reproduce with `python3 -m experiments.v4_run` and `python3 -m unittest discover -s tests -v`. Full records, learned tree parameters, source descriptions, calibration histories and all risk certificates are in [policy_results.json](policy_results.json).


## Follow-up: the planner can now acquire its moments from samples

The [calibration theorem](../../research/v4/calibration_theorem.md) and [four executed cases](calibration_results.json) remove the likelihood table from the learner's interface. It receives sampled target labels and action outcomes, estimates joint moments, and returns a policy with a conditional confidence interval. The two fixed noisy parity worlds, each run with two independent calibration seeds, all yielded exact policy error 9/50 and zero regret. This is one acquired-statistics result; it does not change the null finding for sequential policy distillation.

Each run paid for 576,000 labelled calibration episodes, 576,000 target-label readouts and 921,600 visible readings on conditioned copies. The source's labels, fixed hidden state within an episode, and independent observation noise are substantive assumptions. A standard empirical-moment method ties. The next obligation is to use an acquired representation to reduce further discovery cost.

The expanded suite passes **90 tests**. A separate [persisted-result verifier](../../experiments/v4_verify.py) checked all 448 selected policies, 840 candidate certificates, 32 source likelihood tables and four calibrated policies. [Verification receipt](verification_receipt.json).
