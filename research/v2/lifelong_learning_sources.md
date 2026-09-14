# Lifelong learning and cumulative scientific discovery: bounded source map

This note separates established transfer guarantees from the stronger claim needed for cumulative scientific discovery. None of these papers proves open-ended mechanism discovery, causal abstraction, or universal improvement under arbitrary task sequences.

## Five primary results

1. **Baxter, “A Model of Inductive Bias Learning,” JAIR 12 (2000).** [Primary paper](https://www.cs.cmu.edu/afs/cs/project/jair/pub/volume12/baxter00a.pdf). Baxter models an environment as a distribution over related tasks and a bias learner selecting a hypothesis space from a permissible family. Uniform-convergence Theorem 2 gives task-level and within-task sample requirements in terms of the capacities of the hypothesis-space family and the selected class; a bias that performs well on enough sampled tasks generalizes to new tasks from the *same environment*. The gain is an amortized representation/bias result, not proof that the learner discovers new mechanisms or handles adversarially chosen tasks. The same-environment and finite-capacity assumptions are the essential boundary.

2. **Alquier, Mai, and Pontil, “Regret Bounds for Lifelong Learning,” AISTATS 2017.** [Primary paper](https://proceedings.mlr.press/v54/alquier17a.html). Their EWA-LL meta-algorithm updates a prior over representations as tasks arrive sequentially. Under within-task algorithms having stated regret/statistical guarantees, the meta-algorithm inherits bounds; for a finite predictor family they report an improved `O(1/m)` per-task term (where `m` is samples per task) under their assumptions. Guarantees are expectation-based for general losses and uniform for convex losses, with task data presented in the model's online setting. Transfer is not guaranteed when task structure is unrelated, adversarial, or the representation family omits the useful feature.

3. **Rothfuss et al., “PACOH: Bayes-Optimal Meta-Learning with PAC-Guarantees,” ICML 2021.** [Primary paper](https://proceedings.mlr.press/v139/rothfuss21a.html). Theorem 2 gives a high-probability bound on transfer error over a task distribution in terms of empirical multi-task error, a hyper-posterior/hyper-prior KL term, within-task posterior KL terms, and bounded-loss (or sub-gamma) concentration terms. It is an explicit certificate for meta-generalization under an i.i.d. task distribution and a declared hypothesis/base-learner family. It does not establish sequential causal discovery, and its KL complexity terms make clear that unrestricted growing representations do not receive free generalization.

4. **Amani, Yang, and Cheng, “Provably Efficient Lifelong Reinforcement Learning with Linear Function Representation,” ICLR 2023.** [Primary record](https://www.microsoft.com/en-us/research/publication/provably-efficient-lifelong-reinforcement-learning-with-linear-function-representation/). In a linearly parameterized contextual MDP, tasks may be adaptively chosen from contexts while transition dynamics are context-independent. Under their completeness-style representation assumption, UCBlvd obtains sublinear regret, with reported bound \(\widetilde O(\sqrt{(d^3+d'd)H^4K})\) over `K` task episodes and only `O(dH log K)` planning calls. This is the closest rigorous result to useful learning under adaptive task selection, but the shared linear representation and context-independent dynamics are strong structural assumptions; it is not a theorem for arbitrary physical systems or discovery of new state variables.

5. **Wolpert and Macready, “No Free Lunch Theorems for Optimization,” IEEE TEC 1 (1997).** [Primary DOI](https://doi.org/10.1109/4235.585893). For finite search spaces and the paper's uniform averaging over problem functions (including its time-dependent formulation), algorithms have equal average performance; any advantage on one problem class is offset on another. This does not say useful transfer is impossible. It says that a cumulative-discovery claim must declare the structured task/environment distribution, causal access, or prior; “works on all adversarial physical laws” has no free theorem behind it.

## What is established versus genuinely new

Established: related-task distributions can make learned representations improve future-task generalization; sequential online schemes can bound transfer regret; PAC-Bayes can charge both meta-level and within-task complexity; and shared linear structure can yield sublinear adaptive-task regret.

Genuinely new for the Spectra program: a learner that proposes a physical variable or mechanism, chooses an intervention to discriminate it, produces an independently checkable prediction/certificate, and transfers the mechanism to a new Hamiltonian or coupling while reporting where the abstraction fails. Existing results assume the hypothesis/representation family or task relation sufficiently strongly; they do not prove this discovery loop.

## A valid theorem to pursue here

Fix a finite physical mechanism vocabulary \(\mathcal M\), a bounded intervention set \(\mathcal E\), and a declared distribution over Hamiltonians generated by a sparse latent support rule. Let a learner maintain a representation \(R_t\), choose experiments adaptively, and output a certificate whose test loss is bounded by \(\epsilon\). Prove a decomposition

\[
\text{cumulative test loss}
\leq
\text{from-scratch baseline}
 - \text{transfer gain}(R_t)
 + \text{representation-complexity term}
 + \text{model-misspecification term},
\]

with high probability over independent held-out Hamiltonians and interventions. Require the transfer gain to be positive only after counting discovery experiments, failed hypotheses, certificate verification, and reuse costs. Add an abstention theorem: if no candidate in `M` explains the observations within the declared tolerance, the learner must output a misspecification witness rather than force a mechanism.

This is a finite-family, intervention-aware lifelong-learning theorem. It would be meaningful evidence of compounding scientific capability while remaining compatible with Baxter/PAC-Bayes assumptions and the no-free-lunch boundary. It does not claim universal physical law discovery or adversarial-environment superiority.
