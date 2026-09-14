# Prior art for signed moments, finite-horizon planning, and interaction tests

The signed-moment theorem in `research/v4/signed_moment_theorem.md` is a direct finite-model derivation. The sources below establish adjacent ideas; none is being claimed as a source of that exact formula.

## Three relevant primary sources

1. **Smallwood and Sondik, “The Optimal Control of Partially Observable Markov Processes over a Finite Horizon,” Operations Research 21 (1973), 1071–1088.** [INFORMS primary record](https://doi.org/10.1287/opre.21.5.1071). For finite hidden state, action, observation, and horizon, the optimal value as a function of the belief state is piecewise-linear and convex, represented by finitely many alpha-vectors/conditional plans. This directly supports the theorem’s statement that finite-horizon Bayesian experiment planning is ordinary belief-state dynamic programming. It does not imply that the compressed signed moments are sufficient for arbitrary POMDPs; sufficiency here follows only from the declared binary independent-channel model.

2. **Derezinski, Liang, and Mahoney, “Bayesian Experimental Design Using Regularized Determinantal Point Processes,” AISTATS 2020.** [PMLR primary record](https://proceedings.mlr.press/v108/derezinski20a.html). This paper gives a formal Bayesian experimental-design setting in which a prior over parameters and an information/design criterion define experiment selection, and establishes a connection between Bayesian design and regularized DPPs. It is adjacent prior art for charging experiment-selection cost and optimizing a prior-dependent design objective; it is not a signed-measure or finite-horizon alpha-vector theorem and does not validate the project’s exact risks.

3. **Robin A. A. Ince, “Measuring Multivariate Redundant Information with Pointwise Common Change in Surprisal,” Entropy 19(7), 318 (2017).** [MDPI primary article](https://doi.org/10.3390/e19070318). This work develops an explicit multivariate-information decomposition and discusses XOR/parity-style dependence in which lower-order summaries can fail to capture joint information. The exact two-world construction in the project is simpler and self-contained: two model/sign channels can share all first signed moments while differing in the cross moment `G_AB`, so a one-shot summary cannot determine two-shot planning value. The source supports the general interaction/synergy warning, not the project’s numerical `Q_A` values.

## Scope of the local theorem

For a finite binary experiment family, the signed joint mass of a terminal history is a linear functional of the hidden target sign. Expanding independent raw Bernoulli likelihoods gives multilinear moments of the chosen action probabilities. The displayed two-shot formula and its `G_AB` obstruction are therefore directly derived algebra, not a literature novelty claim. Smallwood–Sondik supplies the general belief-state planning context; the project’s moment representation is a sufficient statistic only under its explicit factorized channel and finite menu assumptions.

The XOR analogy must remain bounded. Univariate marginal summaries can miss a joint interaction, but this does not establish that every representation lacking a particular moment is inadequate: a different sufficient statistic, a joint experiment, or direct posterior computation may recover the same value. The correct conclusion is an obstruction to the specified first-moment summary, not a universal lower bound on planners.

## Valid comparison theorem

For the declared finite model/prior/action menu, compare a moment planner, exact belief-state dynamic programming, and any learned policy under identical likelihood access and shot/action costs. Prove equality of the moment planner and alpha-vector optimum when the required moments are supplied; prove abstention when a withheld cross moment is unavailable or the channel factorization fails. Report planning computation and acquisition costs separately. A tie with the exact Bayesian planner is expected and is evidence of a compact sufficient representation, not algorithmic superiority or open-ended scientific discovery.

Primary records were reopened by the root investigator on 2026-09-10. This corrected an initial agent-provided error in the Ince title, author list, year, volume and DOI. Agent agreement was not accepted as source verification.
