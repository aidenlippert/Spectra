# Exact two-shot signed-moment planning theorem

This note uses the experiment convention in experiments/v4_moments.py. Hidden model m has prior p_m, target sign s_m in {−1,+1}, and binary action a with raw plus probability L_ma in [0,1]:

P(+ | m,a) = L_ma.

Define signed moments

mu0 = sum_m p_m s_m,
mu_a = sum_m p_m s_m L_ma,
G_ab = sum_m p_m s_m L_ma L_mb.

Actions may repeat, and the second action may depend on the first outcome.

## Exact two-shot formula

For first action a, the optimal Bayes risk is

Q_a = 1/2 - 1/2 [
 max_b ( |G_ab| + |mu_a - G_ab| )
 + max_b ( |mu_b - G_ab| + |mu0 - mu_a - mu_b + G_ab| )
].

The first maximum chooses b after a plus outcome; the second chooses b after a minus outcome. The maxima are independent and must be summed. The optimal first action is min_a Q_a. This is exactly the formula used by the finite experiment.

For a fixed branch, the four signed joint masses are obtained by expanding the raw Bernoulli factors. After collecting the common terms, the plus branch contributes the pair G_ab and mu_a-G_ab; the minus branch contributes mu_b-G_ab and mu0-mu_a-mu_b+G_ab. Bayes success is one half plus one half the sum of the relevant absolute signed masses; the displayed Q follows. A zero-probability branch has zero signed mass and permits arbitrary action selection. Repeated actions are included by allowing b=a. An arbitrary prior target bias is represented by mu0.

If instead one uses centered coordinates C_ma = 2L_ma - 1, then the corresponding moments are related by

mu_a^C = 2mu_a - mu0,
G_ab^C = 4G_ab - 2mu_a - 2mu_b + mu0.

The centered expansion has explicit factors of 1/4 in each joint outcome. Mixing centered moments with the raw-probability convention is the source of apparent factor errors; the raw formula above is the convention for this experiment.

## Proof

For outcome x in {+,−}, the signed joint mass is S(x,y) = sum_m p_m s_m P(x|m,a)P(y|m,b). For the plus branch, the two y masses reduce, under the raw convention, to G_ab and mu_a-G_ab. For the minus branch they reduce to mu_b-G_ab and mu0-mu_a-mu_b+G_ab. For any terminal outcome, the Bayes-optimal sign has signed advantage equal to the absolute signed mass. Summing the two absolute masses in each branch and choosing the best b independently gives the stated risk. This also proves sufficiency of (mu0, mu_a, mu_b, G_ab) for exact two-shot planning.

## First moments are insufficient

Take four equally likely models with signs (+,+,−,−). Let raw action A probabilities be (1,0,1,0). Compare:

* World X: raw B probabilities (1,0,1,0), giving mu0=mu_A=mu_B=0 and G_AB=0.
* World Y: raw B probabilities (1,0,0,1), giving the same three first moments but G_AB=1/4.

In World X the formula gives Q_A=1/2. In World Y it gives Q_A=0: the joint two-shot outcome identifies the sign. Thus first signed moments do not determine planning value; the relevant cross moment is necessary. The channels are valid deterministic Bernoulli experiments, including zero-probability outcomes.

## General horizon

For an h-shot sequence, expand the product of raw Bernoulli factors. The signed joint mass is a multilinear polynomial in products of L_ma. Exact finite-horizon planning therefore needs signed moments M_T = sum_m p_m s_m product_{i in T} L_{m a_i} up to order h for the action sequences considered. Dynamic programming over histories is the ordinary finite alpha-vector/POMDP calculation; this moment representation is a sufficient statistic, not an unrestricted discovery theorem.

## Identifiability boundary

The theorem is an exact certificate for the declared finite model family and menu. It does not show that a learner discovers the required cross moments, that the prior or channels are physically correct, or that the representation transfers outside the family. A feature-discovery test must withhold candidate cross moments, prevent hidden model labels from entering fitting, and compare against a Bayesian planner given the same prior and channels. No learner can be required to beat every algorithm: any policy learner can be emulated by an algorithm containing its data and update rule.
