# A bounded constructive proof of complementary scientific learning

**Result.** Within an explicitly supplied family of physical response laws and experiments, one learned law enables acquisition of a second law. Both laws transfer to previously unmeasured controls. Their joint reduction in prediction error is strictly larger than the sum of their separate reductions. A generalization has successive marginal gains that double. Finite noisy acquisition preserves a positive margin, and all acquisition shots are counted.

This is a mathematical construction with a reproducible simulator. It is not a proof of unrestricted scientific intelligence, a new physical law, a quantum computational advantage, or algorithmic superiority. The states are separable and classically simulable. Their value here is that the measurement probabilities have a valid physical realization rather than an invented reward function.

The earlier many-body experiment remains a null result for the stronger intelligence claim. This construction establishes a narrower claim with a precise definition; it does not retrospectively turn ordinary support reuse into a discovery of a research method.

## 1. The capability being proved

Let `R_0, R_A, R_B, R_AB` be downstream error under four information states: neither response law retained, only A retained, only B retained, and both retained. Define complementary value

\[
J=R_A+R_B-R_0-R_{AB}.
\]

A positive value says that acquiring B after A provides more error reduction than acquiring B without A. Equivalently, the joint gain exceeds the sum of the gains in isolation.

This definition alone is insufficient: wrong rules can sometimes produce positive `J`. The construction therefore also requires correct law acquisition with explicit coverage, improved absolute target accuracy, transfer to unused contexts, and a disclosed comparator class. It separately tests whether acquisition costs amortize. No single scalar is used as a certificate of intelligence.

## 2. Physical source and action menu

Let `P_0=X`, `P_1=Z`. For hidden axes `a,b` and target sign `s`, admit the two-qubit source

\[
\rho_{s,a,b}=\frac{I+s v P_a\otimes P_b}{4},
\qquad s\in\{-1,+1\},\quad 0<v\le1.
\]

The nontrivial Pauli product squares to identity and is traceless. The state has eigenvalues `(1+v)/4` and `(1-v)/4`, each twice, so it is positive and normalized. For `v<1` it is the Gibbs state of `H=-sJ P_a tensor P_b` when `tanh(beta J)=v`.

An explicit preparation is a mixture of product eigenstates. Draw a parity bit `Q` with mean `s v`; draw one independent uniform eigenvalue `lambda_1`; set `lambda_2=Q lambda_1`; prepare the two local eigenstates along their source-selected axes. This is a valid preparation in the admitted model. The learner is not given those axes. The source has an unknown stable response to public controls, just as an apparatus can have an unknown calibration law.

The allowed target experiment is one product-Pauli measurement, action `(i,j)` in `{0,1}^2`, retaining its parity outcome `Y`. Pauli orthogonality gives

\[
\Pr(Y=+1\mid s,a,b,i,j)
=\frac{1+s v\mathbf1[i=a]\mathbf1[j=b]}2.
\]

A matched action contains sign information; an unmatched action gives a fair coin. The theorem is relative to this menu. Arbitrary collective measurements, alternative source access and free calibration change the problem.

The resource ledger counts source preparations and measurements as primitive operations. It does not supply a laboratory implementation, energy budget, gate-error model or hardware timing theorem. Consequently this is not resource-complete physical engineering in the sense demanded by the ten endpoint targets.

## 3. Exact optimal one-shot risks

Under a uniform independent prior on the unknown axes and sign, every action that agrees with already known axes has the same probability of matching the remaining axes. Randomizing unknown action axes uniformly gives the same expected risk pointwise for any fixed true axes.

With no law, one law, or both laws retained, the match probability is respectively `1/4`, `1/2`, or `1`. The decision `estimated s=Y` has matched risk `(1-v)/2` and unmatched risk `1/2`. Hence

\[
R_0=\frac12-\frac v8,\qquad
R_A=R_B=\frac12-\frac v4,\qquad
R_{AB}=\frac12-\frac v2.
\]

For a fixed one-shot action, summing the smaller of the two sign likelihood masses for each outcome proves Bayes optimality. Randomizing over actions cannot improve their minimum risk. Actions contradicting known axes give no information and are included in the exhaustive check.

Therefore

\[
J=\frac v8>0.
\]

At `v=1/2`, the exact errors are `7/16`, `3/8`, `3/8`, `1/4`, and `J=1/16`. B saves `1/16` before A and `1/8` after A: its marginal benefit doubles.

The program generates an all-action Bellman certificate for the finite Bayesian decision problem. A separate verifier checks every likelihood, posterior, terminal sign risk and action minimum. It rejects altered likelihoods, omitted actions, duplicate or unreachable nodes and malformed horizons. This is an exact arithmetic check of finite instances, not a formal proof-assistant verification of the general theorem.

## 4. Response laws are acquired from noisy observations

The hidden axes obey unknown binary linear laws

\[
a(x)=m_A\cdot x\pmod2,\qquad
b(z)=m_B\cdot z\pmod2,
\]

where the masks have `d` bits. The linear class is supplied; the masks are not. A rank-`d` calibration design identifies each law, and its coefficients determine predictions on all binary controls. This is structured generalization, not exact-answer lookup. It remains ordinary closed-class identification rather than discovery of an unknown law vocabulary.

Only A has direct calibration access in the implemented protocol: each unit-vector control produces its hidden A bit through independent BSC crossover `eta_A=1/10`.

After learning A, select an A control that was never calibrated, predict its axis, and align the first component of a known-positive two-qubit reference. Measure its second component along axis zero. Conditional on A being correct, the plus probability is `(1+v)/2` if the unknown B bit is zero and `1/2` if it is one.

Postprocess a plus result to bit one with probability `q=v/(2+v)`; a minus result always produces one. Exact arithmetic gives

\[
\Pr(1\mid b=0)=\frac{1-v+q(1+v)}2=\frac1{2+v},
\]

\[
\Pr(1\mid b=1)=\frac{1+q}2=1-\frac1{2+v}.
\]

Thus A compiles a B calibration channel with crossover `eta_B=1/(2+v)`. At the tested visibility, `eta_B=2/5`. The B learner receives these sampled bits, never the simulator mask. No direct B calibration call is used.

This randomized conversion is convenient for a symmetric-channel proof, not acquisition-optimal. A conventional designer can use the same conversion or exploit the original asymmetric observations more efficiently.

## 5. Exact finite-shot coverage

With odd `k` independent observations of a bit through crossover `eta<1/2`, the majority-error probability is

\[
T(k,\eta)=\sum_{i=(k+1)/2}^k {k\choose i}\eta^i(1-\eta)^{k-i}.
\]

A union bound over `d` calibration contexts gives `d T(k,eta)`. For rational noise this is computed with exact integers and fractions. A Hoeffding bound proposes a sufficient shot count; the exact tail independently checks its adequacy. The checker also checks record hashes, rank, majority constraints and inferred coefficients.

Let `G_A` be correct A acquisition. If `Pr(not G_A)<=delta_A` and `Pr(not G_B | G_A)<=delta_B`, then

\[
\Pr(\neg(G_A\cap G_B))\le\delta_A+\delta_B=\delta.
\]

No independence between the two learning stages is assumed. Independence of shots conditional on the preceding correct laws is required. Confidence here is repeated-procedure coverage under the admitted source; it is not posterior certainty that this particular observed model is true.

For `d=6`, the protocol uses 29 shots for each A context and 435 for each B context: 174 single-qubit preparations and 2,610 two-qubit references, totaling 2,784 acquisition shots. The exact summed error bound is approximately `7.458654e-5` per world, smaller than the allocated `0.002`.

## 6. The complementary gain survives acquisition error

At a fixed target, let `c_A,c_B` indicate whether the learned axis prediction is correct. Averaging only over unknown-axis randomization and target outcomes gives

\[
R_A=\frac12-\frac{v c_A}4,\quad
R_B=\frac12-\frac{v c_B}4,\quad
R_{AB}=\frac12-\frac{v c_Ac_B}2.
\]

The one-law risks can never be lower than their ideal values. The both-law risk lies between `(1-v)/2` and `1/2`. It follows that

\[
\mathbb E R_{AB}\le\frac{1-v}2+\frac v2\delta,
\qquad
\mathbb E J\ge\frac v8-\frac v2\delta.
\]

At `v=1/2`, these are approximately `R_AB<=0.250018647` and `J>=0.062481353`.

The exact four-case table is stronger: `J=v/8` when the correctness indicators agree and `J=-v/8` when they differ. It yields `J>=v/8-v delta/4`. The result files deliberately retain the conservative bound above. The same table explains why positive J alone is not an accuracy certificate: both wrong rules also have positive J.

## 7. A chain of discoveries

For an `r`-qubit parity source

\[
\rho=2^{-r}\left(I+s v\bigotimes_{j=1}^r P_{a_j}\right),
\]

retaining `k` correct response laws gives optimal one-shot risk

\[
R_k=\frac12-\frac{v}{2^{r-k+1}}.
\]

Its marginal improvements satisfy

\[
\Delta_k=R_k-R_{k+1}=\frac{v}{2^{r-k+1}},
\qquad \Delta_{k+1}=2\Delta_k.
\]

Separately preparable prefix reference sources let already learned laws align the next calibration, reproducing the same BSC channel. These sources are a substantive access assumption: tracing out sites of the full parity state erases its correlation and does **not** produce the needed prefix reference.

Finite-shot acquisition uses `O(r d log(rd/delta)/v^2)` reference shots for fixed initial-channel quality, and `O(rd)` bits for the learned masks. Prefix sizes make total qubit preparations as large as `O(r^2 d log(rd/delta)/v^2)`. Retaining every raw record costs more than storing masks. The source eigenstate preparation requires linear work in its prefix length in the admitted primitive model.

The increasing marginal result is exact, but the earliest margins are exponentially small. Certifying all of them with imperfect acquisition or source mismatch requires correspondingly smaller errors. This is not an unbounded practical acceleration theorem. The implementation checks the complete one-shot action menu through `r=6`; the written argument covers general finite r.

## 8. Target-only lower bounds and amortization

If hidden axes are uniform across `2^r` possibilities and every unmatched action gives a fair coin, any adaptive target-only policy has probability at most `t/2^r` of hitting the hidden action within `t` shots. Couple its two possible target-sign transcript distributions to the same null fair-coin path until the first matching query. That null path is independent of the hidden axes; a union bound applies without conditioning on previous failures. The two sign transcripts can differ only on a hit, so their total variation is at most `min(1,t/2^r)`. Hence

\[
R_{\rm target-only}(t)\ge\max\left(0,\frac12-\frac{t}{2^{r+1}}\right).
\]

The exact two-qubit Bellman calculation is sharper: the optimal four-shot risk is `103/256`, greater than the chosen error threshold `3/10`. Therefore any **fixed-horizon** reset target-only procedure meeting that threshold needs a budget of at least five shots per task. The acquired two-law procedure meets it with one target shot, including its acquisition failure bound.

For N tasks, its committed shot budget is `2784+N`; the restricted reset comparator needs at least `5N`. The first strict crossing is N=697. At N=8192 the budgets are 10,976 versus at least 40,960.

This is not a lower bound on expected adaptive stopping time. It is not a comparison with an optimal stateful learner that shares information across tasks or accesses calibration. Such a conventional learner can reproduce this strategy. The source preparations also have different sizes; shot counts are not equal energy costs.

## 9. Experiments and adversarial controls

The registered integrated protocol uses three hidden-mask seeds, separate acquisition and evaluation randomness, and 8,192 target trials per world. Each law has six calibration contexts; the 57 nonzero non-unit controls are held out. All 57 squared control pairs are also evaluated with exact probabilities, independently of sampled test outcomes.

All six masks were recovered in the registered run. The exact risks were the theorem values for every world. Frozen and exact-answer lookup coincide because lookup has zero held-out hits. Structured retrieval and the conventional full-feature learner tie the acquired-law strategy. That tie is expected and prevents an algorithm-superiority claim.

The separate representation-revision experiment starts with affine response laws and probes pair interventions. It detects a first omitted interaction, then a second, and transfers the revised quadratic formula to unused higher-weight controls. A cubic term invisible to every calibration query causes held-out disagreement and abstention. Full remeasurement is charged at each changed system; measured warm-shot reduction is zero. This is an additional model-revision result, not evidence for another complementary gain.

Tests cover wrong-prefix calibration, rank deficiency, insufficient samples, corrupted certificates, unseen-context transfer, cubic misspecification and source-model restrictions. Perfectly correlated shot errors and unmodeled source discrepancies remain explicit theorem failure conditions. [Robustness bounds](compositional_transfer_bound.md) state how much mismatch the positive margin can tolerate when uniform discrepancy bounds are actually available.

## 10. What remains open

The construction proves an existence claim about complementary, experimentally acquired knowledge in a finite admitted family. It does not prove that Spectra autonomously discovers the right physical ontology, learns a better research algorithm, or gains an advantage on broad unfamiliar materials systems after including all apparatus costs.

Related-task generalization and sequential representation transfer already have mathematical foundations; these include [Alquier, Mai and Pontil's lifelong regret guarantees](https://proceedings.mlr.press/v54/alquier17a.html), [PACOH's meta-generalization bounds](https://proceedings.mlr.press/v139/rothfuss21a.html), and [Amani, Yang and Cheng's structured lifelong RL result](https://www.microsoft.com/en-us/research/publication/provably-efficient-lifelong-reinforcement-learning-with-linear-function-representation/). No novelty claim is made for the elementary construction here.

The next harder obligation is to discover a transferable representation or experiment-selection rule in the coupled many-body setting, where the useful factorization is not supplied, and retain a gain against matched stateful conventional methods. The first benchmark's rapid operator-space growth shows why this remains difficult. The finite proof here is a reference standard for claims and accounting, not a resolution of that open obligation.
