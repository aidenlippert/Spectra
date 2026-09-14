# Amortized data-learned experimental policies

This note defines a nontrivial finite target for learning experimental strategy. It does not claim unrestricted scientific intelligence. The target is a policy learned from data that transfers across new physical task contexts, amortizes its teacher or planning cost, and preserves a certified risk bound.

## Finite task and experiment model

Let M be a finite set of physical task models with known prior p(m). Let E be a finite menu of available experiments. An experiment e has a finite outcome alphabet Y and a known conditional channel Q_e(y | m). A task context x specifies a posterior over M, a loss for the final decision, a finite horizon H, and a per-experiment cost c(e). Histories are finite strings of chosen experiments and outcomes. A policy pi maps each history and context to either another experiment or a terminal action.

The teacher T is the exact finite-horizon Bayes-optimal policy: dynamic programming over the finite history tree minimizes expected terminal loss plus experiment cost. Its value V_T(x) and risk R_T(x) are therefore exactly computable for this admitted family. The teacher is a reference for planning quality, not an oracle for data outside M or E.

An estimator receives an acquisition data set A of task episodes generated from the declared prior and public channels. It outputs a decision tree or decision list L_A. The output must be data-learned: branch predicates, action choices, and stopping rules are fitted from A, rather than manually inserted from the hidden model labels. A separate certificate checker receives only the candidate policy, public Q, costs, horizon, and a declared promise set P of posterior states or finite parameter cells. It enumerates every state/cell reachable under L_A and checks that the claimed risk bound holds there. Continuous parameters require a finite certified partition and interval error bound; an empirical sample alone is not a certificate.

## Certified transfer claim

For every held-out context x in P, require

R_{L_A}(x) <= R_T(x) + epsilon

and a valid checker certificate with failure probability at most delta only where the channel or parameter uncertainty itself is explicitly included. Held-out contexts must be sampled independently after L_A and its certificate are frozen. No hidden posterior state, teacher action, simulator seed, or evaluation label may enter fitting, branch selection, or stopping.

The transfer metric is total cost, not only terminal accuracy:

C_L(n) = C_acquire + C_fit + C_certificate + sum_{i=1}^n C_execute(L_A, x_i).

Include planning and teacher costs in the matched baseline. If the teacher is run anew for each task, its cost is n C_T(x); if a baseline compiles a fixed policy once, charge its compilation cost once. Report accuracy, certificate coverage, wall-clock or operation cost, and abstentions separately.

## General amortization lemma

Let g be the per-task cost gap between a certified learned policy and a matched baseline on an independently sampled held-out family, with the same risk tolerance. Suppose the learned policy has expected reuse cost at most k_L and the baseline has cost k_B, with g = k_B - k_L > 0. Let D be all acquisition, fitting, certification, storage, and failed-candidate cost paid before reuse. Then

E[C_L(n) - C_B(n)] <= D - n g.

Therefore the learned policy has negative expected excess cost once n > D/g. With a high-probability claim, replace g by a lower confidence bound g_minus and add concentration error for the n held-out task costs. The lemma is elementary accounting; it does not prove that g is positive. Positivity must be established on fresh contexts and at matched certified risk.

## Two-discovery compounding gate

Let A1 produce policy L1 and A2 be acquired after observing failures or residuals on new contexts, producing L2. A2 must add at least one data-selected conditional branch whose predicate is absent from L1 and whose action changes the experiment sequence on a held-out context. Require:

1. Independent certificates for L1 and L2 over the complete finite posterior-state or parameter-cell domain.
2. On a fresh task family E2, L2 preserves the risk tolerance and reduces total cost relative to L1 and all declared baselines.
3. On a second fresh family E3, the branch remains useful under a new coefficient, coupling, or prior context, with acquisition cost for A2 included.
4. The joint gain exceeds the sum of isolated gains under a preregistered value function, or at minimum yields a strictly positive second-stage marginal gain after charging A1.

This is evidence that the learner acquired a reusable conditional operation. It is not evidence of new physical law unless the branch changes the admitted model family or predicts an intervention outcome that was not encoded in the original family.

## Required comparators and separations

* **Exact retrieval cache:** stores solved task-to-policy mappings and returns a match only for repeated contexts. It should fail or abstain on unseen contexts.
* **Nearest-policy retrieval:** retrieves the closest prior context under a public metric and may adapt it, with retrieval and adaptation charged.
* **Stateful Bayesian greedy or optimal:** receives the same prior, channels, and current posterior and recomputes the finite-horizon action online. This is the strongest conventional planning comparator inside the family.
* **Policy compilation/distillation:** receives the same training episodes but is allowed to compress teacher trajectories into a fixed policy. It controls for gains caused only by replacing repeated dynamic programming with a compiled lookup or tree.
* **Scratch learner:** fits from the current acquisition data with no retained policy branches.

Cache reuse proves memory. Distillation proves policy compilation. A data-learned branch that transfers to a genuinely new context proves a stronger research operation. “New physics” requires an expanded model or observable vocabulary and a certified prediction on an intervention withheld from acquisition; a new branch within a fixed supplied menu is not that claim.

## Impossibility boundary

No learner can be required to beat every algorithm or retrieval system. Any learner can be emulated by an algorithm that contains its data, update rule, and policy state. A meaningful claim therefore fixes the finite family, resource model, information available to each comparator, and certificate protocol. The result established here is deliberately conditional: it is a theorem-sized target for amortized strategy learning, not a silent redefinition of unrestricted scientific intelligence.
