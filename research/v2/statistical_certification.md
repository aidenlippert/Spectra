# Finite-shot certification for sequential mechanism discrimination

This note gives a small, explicit certificate suitable for adaptive experiments. It separates Bayesian posterior scores from frequentist anytime-valid evidence. The latter is needed when the experimenter may stop, peek, switch mechanisms, or reuse the same protocol many times.

## Bernoulli hypotheses

For one fixed experiment, observe conditionally independent outcomes (Y_t\in\{0,1\}). Under simple model (M_j), (Y_t\sim\mathrm{Bernoulli}(p_j)). With prior weights πj, the finite-hypothesis Bayesian posterior is

\[
\Pr(M_j\mid Y_{1:n})=\frac{\pi_j p_j^{S_n}(1-p_j)^{n-S_n}}{\sum_k\pi_k p_k^{S_n}(1-p_k)^{n-S_n}},\qquad S_n=\sum_tY_t.
\]

This is useful for ranking hypotheses, but a posterior threshold is not a distribution-free error certificate: it depends on the prior and on the hypothesis class containing the truth. A misspecified class can assign high posterior to the least-wrong model.

For a frequentist sequential certificate, test each null (M_j) against a fixed alternative (M_k) with the likelihood ratio

\[
L_n^{k:j}=\prod_{t=1}^n\frac{p_k^{Y_t}(1-p_k)^{1-Y_t}}{p_j^{Y_t}(1-p_j)^{1-Y_t}}.
\]

Under (M_j), (L_n^{k:j}) is a nonnegative martingale with mean one. Ville's inequality gives
\[
\Pr_{M_j}\!\left(\sup_nL_n^{k:j}\ge1/\alpha\right)\le\alpha,
\]
so stopping at any data-dependent time after crossing (1/\alpha) controls the type-I error. This is the e-process principle developed in [Howard et al.](https://arxiv.org/abs/1808.03204). If the alternative is composite, average likelihood ratios over a predeclared alternative prior; the resulting mixture remains a nonnegative supermartingale under the null.

For K candidate mechanisms, allocate αj with Σjαj≤α and reject candidate j only when its null e-process crosses (1/αj). A union bound controls the probability of any false rejection by α. Adaptive experiment selection is valid when the selected experiment and its next betting factor are predictable from past data and the null specifies the conditional outcome law.

## Explicit sample costs

If two Bernoulli means differ by Δ=|pk−pj|, Hoeffding gives
\[
\Pr(|\hat p_n-p|\ge r)\le2e^{-2nr^2}.
\]
Thus a fixed-α interval of half-width r needs (n\ge\log(2/\alpha)/(2r^2)). To identify the larger mean with gap Δ and error at most α, use (r=Δ/2):
\[
n\ge\frac{2}{\Delta^2}\log\frac{4}{\alpha}.
\]
For K mechanisms and a simultaneous certificate, replace α by α/K (or use the exact likelihood-ratio allocation). This is conservative but explicit and remains valid under optional stopping when replaced by a confidence sequence, such as those constructed by [Waudby-Smith and Ramdas](https://arxiv.org/abs/2010.09686).

Likelihood ratios give the sharper scale. Under the true (M_k),
\[
\mathbb E_k[\log L_n^{k:j}]=nD(\mathrm{Bern}(p_k)\Vert\mathrm{Bern}(p_j)).
\]
Consequently the expected crossing time is approximately ​(log(1/\alpha_j)/D(p_k\Vert p_j)), up to overshoot and transient terms. No finite-shot method can guarantee uniform fast discrimination when the KL divergence tends to zero. For small gaps away from 0 and 1, (D(p_k\Vert p_j)\approx\Delta^2/[2p(1-p)]), recovering the inverse-square scaling.

## Version-space elimination and e-process use

Maintain the candidates whose null e-processes have not crossed their rejection thresholds. At each allowed stop, report the surviving version space, the crossed evidence, and the assumptions governing each likelihood. An experiment is informative when it maximizes a predeclared lower bound on the minimum pairwise KL divergence among survivors; the choice may use past data but not the next outcome. If all candidates are rejected, report model-class failure and select a new mechanism language. Do not convert “best surviving posterior” into a discovery certificate.

For a quantum simulator, each candidate predicts a Bernoulli event after a fixed POVM and control sequence. Calibration data must be either included in the likelihood or held out. If the simulator prediction has total-variation error at most ε uniformly for that experiment, then the effective observable gap can shrink from Δ to at least Δ−2ε. A Hoeffding guarantee therefore uses Δeff=max(0,Δ−2ε), and a likelihood certificate must use a robust null set containing every (p_j'∈[p_j−ε,p_j+ε]). If Δ≤2ε, that experiment cannot certify discrimination regardless of shot count. Unknown calibration error must be estimated with its own confidence sequence or incorporated as a nuisance parameter; otherwise the statistical certificate only certifies the simulator, not the physical mechanism.

The resulting artifact is an anytime-valid finite-shot diagnostic, not proof that a mechanism is physically complete. Its strongest negative output is either a rejected candidate under controlled error or an explicit failure of the admitted model class.

## Noisy transfer and calibration accounting

For the proposed transfer experiment, keep the ideal Bayesian table separate from the calibration theorem. Let source A be a direct parity source passed through a BSC with error eta_A=1/10. Let the learned B channel use the aligned two-qubit reference with s=+1, j=0, and source visibility v=1/2. The stated post-processing gives a BSC with error
\[
\eta_B=\frac{1}{2+v}=\frac25,
\]
because the plus branch maps to bit 1 with probability (r=v/(2+v)=1/5), while the minus branch maps to bit 1 deterministically. These are conditional channel parameters; they should not be silently treated as unconditional frequencies after selecting calibration successes.

For each of the (d=6) recorded masks, use a fresh calibration sample and an anytime-valid confidence sequence for the relevant Bernoulli parameter. Allocate δ_A=δ_B=0.001 per record, or use a single familywise allocation with total δ≤0.01. A union bound then gives an event E, with probability at least 1−δ, on which every acquired A and B map obeys its declared error bound. Conditional on E, the exact ideal Bayesian experiment tables may be applied to fresh evaluation worlds, provided the evaluation masks and contexts are sampled from the declared prior independently of the calibration outcomes.

This conditioning requirement matters. Selecting or retaining masks because their calibration confidence intervals look favorable can change the distribution of unretained masks and invalidate an unconditional fixed-prior posterior calculation. The safe alternatives are to freeze the map after calibration and evaluate on independent fresh worlds, include calibration uncertainty as a nuisance parameter, or average analytically over the finite prior. For the small GF(2) family, exact enumeration over all masks and hidden-axis configurations is preferable to treating a Monte Carlo table as an exact theorem. Monte Carlo can still estimate a fixed-world risk, but its confidence interval must be clustered by world/map and must include calibration randomness.

If the ideal complementarity score is (J_\mathrm{ideal}) and the score is formed from four bounded Bayes-risk terms whose total range is at most one per term, the conservative failure-event perturbation is
\[
J_\mathrm{actual}\ge J_\mathrm{ideal}-4\delta
\]
in expectation (or with the corresponding high-probability additive term if each risk estimate is separately confidence-bounded). With (v=1/2), any claimed positive lower bound such as (v/8) must therefore retain an explicit calibration subtraction; δ=0.01 leaves a positive bound only when the ideal margin exceeds 0.04. This is a robust transfer statement, not evidence that the learned map is correct outside the declared mask/context prior.
