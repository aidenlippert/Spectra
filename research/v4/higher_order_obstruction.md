# A hierarchy of missing joint information

This is a directly derived finite statistical obstruction. No novelty claim is made. It complements the **two-shot, three-hidden-bit** construction in `tests/test_v4_feature_obstruction.py`; that construction is not a three-shot test.

## Two-shot action-selection obstruction

Let u,v,w be independent uniform bits and the target sign be s=(-1)^(u+v). Four deterministic binary actions have outputs:

| World | A | B | C | D |
|---|---|---|---|---|
| Left | u | v | w | w |
| Right | w | w | u | v |

Use the same uniform eight-model prior and equal preparation metadata in both worlds. Each individual action has the same joint distribution with the target: an independent fair bit. It also produces the same posterior-mass multiset over the eight models. Consequently **all twelve implemented per-action features are identical**, including sign risk, model information gain, likelihood variance, posterior purity, and preparation weight. The test checks equality of the entire feature matrices.

In Left, beginning with A or B allows the other member of that pair to reveal the target parity; beginning with C or D leaves only one useful bit available in the last shot. Thus the exact first-action risk vector is (0,0,1/2,1/2). In Right it is (1/2,1/2,0,0).

Any rule that sees only the feature matrix, including any randomized rule with the same training history, chooses the same distribution over first actions in both worlds. Its mean regret under a uniform draw between the worlds is exactly 1/4. At least one world's regret is therefore at least 1/4. The exact planner, which sees the cross-action alignment of likelihoods, achieves zero in both.

This is an information obstruction for the specified feature interface. It does not show that these deterministic channels belong to the narrower noisy Clifford benchmark, nor that this obstruction alone caused its empirical null. More training or a larger regressor cannot repair this interface uniformly over the larger finite-channel class. A joint-action feature, retained likelihood alignment, or an additional planning query can repair it.

## General horizon theorem

For h>=2, draw hidden bits x_1,...,x_h,z uniformly **once per prediction episode**. Action i reveals x_i. All actions within that episode refer to the same hidden tuple. Compare world P, with sign s=(-1)^(sum x_i), and world I, with sign s=(-1)^z.

For every proper subset S of the h visible bits,

E_P[s product_{i in S} (-1)^x_i] = E[product_{i not in S} (-1)^x_i] = 0.

The corresponding signed moment in I is zero by independence of z. Expanding raw probabilities L_i=x_i in centered variables shows that all raw signed likelihood moments of degree below h also agree. Repeated indices cannot evade the obstruction: deterministic L_i satisfies L_i^k=L_i, so repetitions reveal no new bit. All unsigned moments, at every order, agree because the observation law is identical in both worlds.

After h distinct actions, P has risk zero, whereas I retains risk 1/2. A complete order-h signed moment representation can distinguish these planning values. Moments of order below h cannot. This proves necessity of information of arbitrarily high interaction order for a family with growing horizon; it is not a lower bound on every possible representation or algorithm.

The h=2,...,5 cases are also checked by exact finite enumeration in the test file, independently grouping terminal outcomes to compute Bayes risk.

## Physical source and resource boundary

A commuting realization stores the tuple in computational-basis qubits and measures each visible bit in Z. The target register is inaccessible during prediction. A fresh tuple is drawn between episodes. Within an episode one may nondestructively read different commuting registers, or prepare fresh copies conditioned on the same fixed tuple. **Drawing fresh independent hidden bits for every action would invalidate the parity-recovery conclusion.**

If a labelled calibration apparatus additionally supplies s and all h action outcomes per independent episode, the statistic Z=s product_i (-1)^x_i lies in [-1,1]. Its mean is 1 in P and 0 in I. For N independent episodes, Hoeffding gives Pr(|mean(Z)-E Z|>=epsilon)<=2 exp(-N epsilon^2/2). Thus N>=2 epsilon^(-2) log(2/delta) suffices to estimate this moment, charging h visible measurements and the target-label preparation/readout per episode. This source of labels is an extra declared resource, not free discovery access.

The theorem neither proves efficient acquisition of a general order-h tensor nor demonstrates a learner discovering that tensor. Model misspecification, unavailable labels, and growing action menus remain separate obstacles. See [primary-source context](moment_sources.md).
