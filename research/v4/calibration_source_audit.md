# Calibration-source audit for signed-moment estimation

This protocol removes the likelihood table from the learner. The learner receives labelled episodes consisting of a calibrated target sign and sampled outcome bits, then estimates the moments used by the signed-moment planner. The audit is about statistical validity; it does not establish scientific discovery or superiority over a conventional estimator.

## Correct episode structure

Let hidden model m be drawn once from the task prior p. The calibrated preparation supplies its target label s in {−1,+1}. For each chosen action a, obtain an outcome bit X_a in {0,1} from an independent copy conditional on the same hidden m. Conditional independence of the copies gives

E[s] = mu0,
E[s X_a] = mu_a = sum_m p_m s_m L_ma,
E[s X_a X_b] = G_ab = sum_m p_m s_m L_ma L_mb.

The diagonal term G_aa must use two independent copies with the same hidden m. Reusing one Bernoulli outcome would give E[s X_a^2] = E[s X_a] = mu_a, not G_aa. A labelled episode may therefore contain a source preparation, two or more independently measured copies, and the action choices applied to those copies. The source label, source-preparation cost, copy-preparation cost, measurement shots, and failed calibration attempts are all resource costs.

The crucial distinction is shared hidden identity. If a fresh m is redrawn independently for the a and b measurements, the product estimates a product of marginal quantities, not G_ab. For example, it estimates E[s X_a] E[X_b] if the label is attached only to the first draw, or another model-dependent product if labels are attached separately. Neither equals E_m[s_m L_ma L_mb] in general. A mixture-redrawn experiment can therefore pass ordinary sampling checks while invalidating the planner's cross moments.

The protocol must state whether the source creates one persistent hidden model for the episode, whether copies are conditionally independent given that model, and whether an action changes the hidden model. A stateful source or action-dependent transition requires a different moment theorem.

## Finite-sample bound

For A actions, estimate the K = 1 + A + A(A+1)/2 quantities consisting of mu0, all mu_a, and symmetric G_ab with a <= b. Use N independent labelled episodes for each quantity, with fresh source/copy randomness. Every summand lies in [-1,1]: use s for mu0, s X_a for mu_a, and s X_a X_b for G_ab. Hoeffding's inequality gives

Pr(|hat(theta)-theta| >= epsilon) <= 2 exp(-N epsilon^2 / 2).

By a union bound, all K estimates are within epsilon with probability at least 1 - delta whenever

epsilon >= sqrt((2/N) log(2K/delta)).

If different moments use different sample counts N_q, replace the common epsilon by the corresponding per-moment radii and union the individual failure probabilities. If the action menu or moment set is selected after seeing data, either hold out a fresh calibration set or use a simultaneous bound over the full predeclared menu; post-selection without such a correction invalidates the stated confidence.

The planner should propagate these intervals through its rational risk formula, using interval arithmetic or a robust worst-case optimization over the moment box. It must abstain when action rankings overlap within the certified uncertainty or when model/label uncertainty is not included. A point estimate followed by an unqualified action choice is not a certificate.

## What this audit establishes

The source protocol establishes an unbiased finite-sample route to the raw signed moments without exposing likelihood parameters L_ma or hidden model identities to the learner. It can support an exact or interval-certified two-shot planner under the admitted finite model and source assumptions.

It does not establish two successive scientific discoveries. The action menu, target-sign semantics, source distribution, and moment feature family are supplied. Discovering a useful branch within this menu is policy learning or estimator reuse unless the learner also acquires a new representation or intervention rule and transfers it to an independently generated family.

It also does not establish an advantage over standard statistics. A conventional Bayesian or method-of-moments estimator given the same labelled episodes, prior, action menu, and confidence budget is a matched comparator. Any claim of lower cost must include source-label acquisition, copy preparation, measurement, fitting, certification, and failed or abstained policies. The learner must be evaluated on held-out source episodes or new task families, not on the same labelled records used to choose moments.

## Necessary assumptions and failure tests

The guarantee requires iid episodes from the declared prior, accurate and independent source labels, stationarity across calibration and evaluation, finite known action menu, bounded outcomes, and conditional independence of copies given m. It also requires that the source label not leak hidden model metadata unavailable at deployment. Label noise must be measured and included as a channel; an unknown systematic label bias can invalidate every signed moment. Drift, common-mode detector noise, copy correlation, action back-action, and adaptive stopping need explicit models or holdout correction.

Useful negative controls are: deliberately redraw m between paired copies, permute source labels, replace paired copies by independent mixture draws, and introduce a known label-noise rate. The estimator should reject or widen its interval under these violations. A successful numerical fit under the wrong redrawn-mixture protocol is evidence of implementation consistency only, not evidence for G_ab or for reusable scientific structure.
