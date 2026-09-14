# Adaptive experimental policy: bounded transfer test

This benchmark tests transfer of an acquired research rule using a declared
finite quantum source. It does not claim discovery outside that model class.

For (a,b\in\{X,Z\}), sign (s\in\{\pm1\}), and visibility (v=1/2), the
source state is

\[
\rho_{s,a,b}=(I+s v P_a\otimes P_b)/4.
\]

A (P_i\otimes P_j) measurement returns (+1) with probability

\[
\Pr(+|i,j)=\frac{1+s v\mathbf 1[i=a]\mathbf 1[j=b]}2.
\]

The implementation checks this against an independent (4\times4) density
matrix and Born-rule calculation and samples exact rational probabilities. Each
episode also has a hidden binary calibration map: a binary knob vector is mapped
to a source bit by a private mask and passed through a known bounded binary
symmetric noise channel. Masks are never returned to the learner.

The benchmark has three sequential episodes. Stage 1 learns a public source
relation from calibration and target shots. Stage 2 changes the hidden masks and
holds out probe actions. Stage 3 changes source contexts and coefficients and
uses an independent probe bank. Transfer requires reaching a predeclared
posterior-odds threshold, predicting held-out target probabilities, and using
fewer shots than matched baselines. A misspecified episode adds a source relation
outside the hypothesis grid; correct behavior is abstention or an explicit
request for an expanded model class.

For a finite hypothesis grid Θ over source mechanisms and calibration maps,
maintain (p_t(\theta)). For each admissible action (u), select the action
maximizing expected information gain per cost:

\[
u_t=\arg\max_u\frac{H[p_t]-
\mathbb E_{y\sim p_t(y|u)}H[p_{t+1}(\theta|y,u)]}{c(u)}.
\]

Also compute exhaustive finite-grid dynamic programming for the remaining shot
budget. This is a conditional optimum for the declared grid, not universal
optimality. The posterior and certificate must range over joint
(mechanism, calibration) hypotheses; marginalizing calibration can make two
mechanisms observationally equivalent.

If an action separates two hypotheses by Bernoulli total variation Δ, a
Hoeffding test has error at most {exp(-2nΔ^2)} after (n) independent
shots. If every allowed action gives identical distributions for a pair, no
policy can identify it. Adaptive action-selection computation counts toward the
budget.

Compare fixed probe, random, Bayesian-updating-without-policy-learning,
retrieval-only, transferred-policy, and exhaustive-planning baselines. Match
shot, action, simulator, precision, wall-clock, and policy-computation budgets.
Report shots to threshold, held-out log loss/Brier score, calibration error,
abstention, and training/reuse cost separately. If retrieval matches transfer,
report no compounding gain.

The design is informed by Evans, Harper & Flammia, [Scalable Bayesian Hamiltonian
Learning](https://arxiv.org/abs/1912.07636), Zhang et al., [Identifiability
Guarantees for Causal Disentanglement from Soft Interventions](https://arxiv.org/abs/2307.06250),
and Polyanskiy & Wu, [Dualizing Le Cam's Method](https://arxiv.org/abs/1902.05616).
