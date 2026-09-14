# v4 gate: data-learned experiment selection after sequential discoveries

## Claim boundary

The target is a finite-domain learned policy: after acquiring two physical mechanisms from finite-shot data, the system selects the next experiment in a coupled quantum family. The claim is that a learned conditional rule transfers to withheld contexts and lowers expected downstream risk or experiment cost after acquisition and planning costs are charged. This is a benchmark theorem/measurement, not a claim of general scientific intelligence.

An exact Bayesian planner with the same model, observations, prior, and action menu may tie the learned policy. The benchmark must therefore test transfer and amortized reuse against explicit finite comparators; it must not require beating every conventional algorithm.

## Minimally credible finite benchmark

Use a small coupled Pauli family, for example two or three qubits with a hidden Hamiltonian selected from a declared finite grammar containing local terms, one coupling term, and a continuous bounded nuisance parameter (coupling strength, readout bias, or drift). Generate finite-shot circuit outcomes from an explicit state-preparation, evolution, and measurement model. The learner sees observations and controls, never hidden labels or exact Hamiltonian parameters.

Two mechanisms are acquired sequentially from disjoint calibration episodes. Their outputs must be reusable predicates or sufficient statistics whose validity region is recorded. The downstream task presents new continuous parameter values, new initial states, and at least one new coupling context. These held-out contexts must make exact answer-table lookup impossible.

The teacher is an exact finite-horizon Bayesian planner or dynamic program over the declared finite latent model. Teacher episodes are generated before policy fitting, and teacher planning time, model enumeration, certificate checking, and data generation are charged. The learner receives state/action labels only as training targets, never hidden truth. Train a small decision tree or finite grammar policy on teacher state summaries, with a fixed complexity budget and an abstain action when the summary is outside its certified domain.

## Required comparators

Evaluate the same held-out episodes with:

1. **Frozen:** a stage-one policy and mechanism state, never updated.
2. **Exact retrieval:** lookup of prior state/action pairs; it must miss on withheld continuous contexts.
3. **Stateful greedy:** reuses learned mechanisms but chooses the next experiment by a fixed local information or uncertainty score.
4. **Stateful Bayesian planner:** recomputes the finite-horizon optimal action from the same posterior and action menu.
5. **Learned policy:** executes the distilled rule, including abstention outside its certified domain.

The Bayesian planner is the principal conventional comparator. A tie is expected and does not invalidate the finite transfer claim. Any claimed advantage must be against a named implementation under matched model access, not against an unspecified class of conventional methods.

## Genuine learning gate

A policy counts as data-learned only if all of these hold:

- Its action rule is fit from teacher-labelled episodes and held-out continuous/context variables, rather than supplied as an algebraic formula or selected by a hidden mechanism name.
- A frozen answer cache receives the same serialization and lookup budget and has zero or predeclared low hit rate on evaluation contexts.
- The policy must select actions correctly on contexts whose exact state vector and parameter values were absent from training.
- A feature-erasure ablation removes the learned mechanism summary while retaining raw observations; performance must change according to the preregistered prediction.
- A shuffled-label or constant-policy control fails the task at the predicted rate.
- The rule is evaluated after adding a second discovery; a policy trained only before that discovery is a genuine frozen baseline.

Storing a discovered mask, support, or sufficient statistic is reusable structure and may be valuable, but is not by itself learned experiment selection. A structured compiler that reconstructs the same rule from the same data is an allowed tie. The report must distinguish cache/retrieval, supervised policy distillation, and planner recomputation.

## Predeclared positive result

Before seeing evaluation outcomes, specify a risk tolerance, confidence level, policy-complexity cap, minimum number of independent seeds, and amortization horizon (N). A positive finite-domain result requires:

1. The learned policy passes an independent action/risk certificate on its declared domain.
2. It transfers to withheld continuous noise/context parameters and new coupled systems within the declared grammar.
3. It beats frozen and exact retrieval at equal total cost, including teacher acquisition, policy fitting, verification, and execution.
4. The comparison with the stateful Bayesian planner is reported honestly; a tie means no algorithmic superiority.
5. The learned policy’s cumulative cost or risk improvement survives the four-way ablation and confidence interval.
6. Out-of-domain inputs trigger abstention or a measured failure, rather than silent extrapolation.

A null result is valid if planner recomputation ties, acquisition does not amortize, the rule fails continuous transfer, or the learned policy is no better than a cache after charging fitting and verification.

## Certificate and leakage requirements

Use independent source seeds for teacher, training, validation, and final evaluation. Hidden Hamiltonians, nuisance parameters, exact optimal actions, and withheld outcomes must not enter feature construction, stopping rules, or certificate thresholds. Teacher labels may train the policy but cannot be used to certify it on the same episodes.

The independent verifier must recompute finite-horizon risk over the declared latent model, integrate or bound the continuous nuisance parameter, and check every action branch of the deployed policy. Monte Carlo estimates may supplement but cannot replace a finite exact or concentration-bound certificate. If a policy is certified only on a finite grid, state that its guarantee is grid-conditional and require an explicit interpolation margin for continuous transfer.

Report acquisition shots, teacher planning operations, learner fitting operations, verifier operations, policy storage, lookup, execution shots, failures, and cumulative cost separately. Do not call equal physical-shot counts equal total cost.

## Interpretation

The strongest defensible outcome is: “A finite data-learned experiment-selection rule, trained from charged teacher episodes, transfers to withheld contexts under an independently checked model-conditional risk bound and amortizes against frozen/retrieval baselines.” It remains compatible with a conventional Bayesian planner tying the result. It does not establish unrestricted causal abstraction, new physical-law discovery, or engineering superintelligence.
