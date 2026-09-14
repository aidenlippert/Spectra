# v4 policy-run adversarial audit

Audited `v4_run.py`, `v4_physics.py`, `v4_planner.py`, `v4_moments.py`, and `v4_learner.py`, then ran `run(small=True)`. The pilot is useful as an exploratory policy-distillation benchmark, but it does not yet establish sequential learned-policy compounding.

## Concrete blockers

1. **The physics and sequential mechanisms are supplied.** `make_family` supplies the finite Clifford gate grammar, candidate pool, likelihood table, model signs, and visibility construction. `context` performs two calibration observations, but the “first” and “second” stages are datasets for policy fitting rather than discoveries acquired by the learner and then used to alter a downstream experiment. The learned object is a regression rule over planner-derived features, not a discovered physical mechanism. The code’s own interpretation correctly says physics and Bayesian likelihoods are supplied.

2. **Teacher targets and features expose the planner’s answer structure.** `features` receives the exact posterior, full likelihood table, and exact action-conditioned statistics. `y` is directly the signed-gram two-shot action value minus one feature. This is legitimate supervised policy distillation, but it must be described as such. It cannot be counted as discovering the response law or experiment design from raw physical data.

3. **The signed-Gram comparator is the dominant exact baseline.** `signed_gram_plan` computes the optimal first action from exact finite moments without posterior child searches. In the small pilot its mean exact regret is zero while learned policies have positive regret; the cumulative, scratch, batch, nearest, and exact-retrieval rules tie each other and can be worse than the feature-erasure control. This means any “learning improves planning” claim currently fails. Root’s removal of redundant child-search accounting is necessary, but the exact signed-Gram baseline must remain primary and its construction cost must be charged.

4. **The apparent cumulative policy gain is absent.** In the small run, `cumulative_learned`, `scratch_second_stage`, `batch_conventional_learner`, `nearest_record_retrieval`, and `exact_retrieval` have identical mean risk/regret. `feature_erasure` is slightly better than the learned policy, showing that the current feature set/model is not evidence of useful sequential complementarity. Report this as a null pilot, not a positive gate.

5. **Costs are not commensurate.** The learned methods include feature extraction and measured prediction/certificate time, while the signed-Gram method includes teacher construction but still rebuilds a policy certificate. The fixed target protocol, teacher episodes, raw calibration shots, model-table generation, and any mechanism-acquisition costs need one common amortized ledger. CPU wall time across Python paths is not a scientific resource theorem without operation counts and identical verification obligations.

## Smallest honest next experiment

Keep the supplied finite physics family and state the narrow claim as policy distillation only. Acquire a mechanism from finite-shot records in stage 1, hide one response parameter or coupling from the policy feature construction, and require stage 2 to choose an experiment whose likelihood cannot be computed from the supplied posterior table without using the retained mechanism. Compare against signed-Gram exact planning, a stateful Bayesian planner, exact retrieval, frozen, and feature-erasure controls. Charge teacher generation, mechanism acquisition, fitting, verification, and execution over a declared target count.

The minimum positive gate is a preregistered reduction in total risk/cost against frozen and exact retrieval on held-out continuous contexts, with no improvement claim against signed-Gram/Bayesian planners unless it is actually measured. If the exact baseline ties, the strongest result is a verified amortized policy representation, not scientific mechanism discovery.

No hidden-truth access bug was found in the evaluator path, but the model family, likelihoods, and action grammar are intentionally supplied. The current run must therefore remain labelled exploratory supervised planning, with no sequential compounding claim.
