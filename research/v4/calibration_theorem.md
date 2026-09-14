# From sampled calibration records to a certified experiment policy

The planner can acquire its signed moments without receiving a hidden-model likelihood table. This construction removes one supplied-information dependency in V4. It is a finite acquisition result, not a demonstration of two successive scientific discoveries or a superior learning algorithm.

## Source and estimator

There are A binary actions and a target sign s. One calibration episode draws a hidden model m from the target prior, supplies a calibrated target label s_m, and performs zero, one or two specified actions on copies conditional on **that same m**. Their measurement randomness is independent given m. Different episodes are independent. The learner receives only labels and binary outcomes. Calibration and prediction share the same distribution and preparation/readout conditions.

For raw plus likelihood L_ma, the random variables s, sY_a and sY_aY_b estimate respectively mu0, mu_a and G_ab. Repeated actions use independent noise on copies of the same m; squaring one observed bit estimates a different object. Redrawing m for each action also estimates a different object. The negative-control test distinguishes a correct diagonal moment 0.41 from the incorrectly redrawn value 0.25.

There are K=1+A+A(A+1)/2 moments. Each sampled statistic lies in [-1,1]. For N independent episodes per moment, Hoeffding and a union bound give

Pr(any moment error > epsilon) <= 2K exp(-N epsilon^2/2).

No independence between different moment estimates is needed for the union bound. The implementation uses an integer t with 2K 2^(-t)<=delta and N=ceil(2t/epsilon^2). Since exp(-t)<=2^(-t), this conservatively certifies the target coverage without treating a rounded floating logarithm as a proof step.

## Uniform policy-risk and regret theorem

Allow error radii e0 for mu0, e1 for every mu_a, and e2 for every G_ab. For a fixed first action a and branch-dependent second actions b_-,b_+, the four signed terminal masses are

(mu0-mu_a-mu_b-+G_ab-, mu_b--G_ab-, mu_a-G_ab+, G_ab+).

The sum of their absolute estimation errors is at most e0+4e1+4e2. For fixed terminal decisions d_h in {-1,1}, true risk is

R(pi)=1/2 - (1/2) sum_h d_h S_h.

Consequently every fixed two-shot policy obeys

|R(pi)-Rhat(pi)| <= kappa := e0/2+2e1+2e2.

The bound is simultaneous over **all** policies once the moment-error event holds, including policies chosen after observing calibration data. No additional union bound over policies is needed.

The signed-moment algorithm exactly minimizes the estimated signed-score objective over first actions, branch-dependent second actions and terminal decisions. If pihat is that minimizer and pi* the true optimum,

R(pihat) <= Rhat(pihat)+kappa <= Rhat(pi*)+kappa <= R(pi*)+2kappa.

With a common radius epsilon, the conditional selected-policy interval is Rhat±(9/2)epsilon and regret is at most 9epsilon. Empirical moments need not correspond to any realizable probability distribution for this proof: the surrogate is a linear signed score. Its value must not be called an exact probability. Intersect the confidence interval with [0,1], and refuse if it is empty.

The test suite independently enumerates all 128 action/decision trees for each of eight small likelihood worlds. It checks the uniform bound, the surrogate minimization and the resulting regret inequality. A second implementation integrates actual joint outcomes to evaluate the learned policy.

## Executed source cases

The source stores uniform hidden bits u,v,w, with target s=(-1)^(u+v). Its actions read (u,v,w,w) or (w,w,u,v), with independent bit-readout visibility v=4/5. The first two useful observations together reveal noisy parity with optimum error (1-v^2)/2=9/50. An irrelevant first observation leaves error 1/2.

A feature-only first-action rule seeing identical individual-action summaries in both worlds has mean error (9/50+1/2)/2=17/50 over an equal world mixture. A learned joint-moment policy can distinguish the useful pair. In four runs—two calibration seeds per fixed world—it selected an optimal pair and achieved exact conditional error 9/50. All observed moment errors and risk intervals satisfied the stated bounds. These four executions do not themselves validate the confidence theorem; that comes from the independent-episode and source assumptions above.

For A=4, epsilon=1/40 and delta=1/100, K=15, t=12 and N=38,400. Each run consumed 576,000 simulated labelled episodes, 576,000 target-label readouts and 921,600 visible readouts on conditionally prepared copies. The conservative union-failure upper bound is 15/2048, below 1/100. The guaranteed regret bound 9/40 is loose; the observed exact regret is zero. The source construction and its hardware implementation cost remain admitted assumptions, and no hardware experiment was performed.

## Exact remaining boundary

This is a standard empirical-moment estimator followed by a programmed sufficient-statistic planner. An equally equipped conventional estimator ties it. The representation is acquired numerically from raw outcomes, but the choice of moment vocabulary and error theorem is supplied by the research implementation. It does not show the system inventing that vocabulary or using one discovery to make a second discovery cheaper.

The result advances the pipeline from supplied likelihoods to sampled calibration with explicit uncertainty. The remaining compounding obligation is positive acquisition and transfer at a further stage, with the same calibrated-source and resource accounting. See [results](../../results/v4/calibration_results.json), [independent bound audit](estimated_moment_bound_review.md), and [source audit](calibration_source_audit.md).
