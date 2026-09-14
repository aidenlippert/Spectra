# Robust compositional transfer

These are model-conditional bounds, not evidence that a real apparatus meets their premises.

## Calibration under bounded discrepancy

Let the ideal processed calibration label be a bit through BSC crossover `eta < 1/2`. Suppose the actual, independently prepared shot distribution is within total variation `epsilon_j` of that ideal distribution for every admitted context and hidden bit. Its error probability is then at most `eta + epsilon_j`. Independence need not imply identical distributions: coupling each error indicator to an independent Bernoulli variable with probability `eta_* = eta + max_j epsilon_j` stochastically dominates the majority-error event.

For odd `k`, the probability of any error among `d` majority labels is bounded by

\[
\delta\le d\sum_{i=(k+1)/2}^k {k\choose i}\eta_*^i(1-\eta_*)^{k-i}
\le d\exp[-2k(1/2-\eta_*)^2].
\]

The bound requires `eta_* < 1/2`. At the compiled B channel with visibility `v`, `eta = 1/(2+v)` and the remaining noise margin is `v/[2(2+v)]`. Source mismatch must be smaller than this margin. Exact binomial arithmetic applies when `eta_*` is rational and `k` is an integer. Shot cost diverges quadratically as the effective margin approaches zero.

Stochastic postprocessing contracts total variation, so a bound on the raw reference distribution is also a valid bound on its processed label. It does not certify the raw source itself.

## Composition through interfaces

Let `K_j` and `Ktilde_j` act on the same declared interface state spaces, with a uniform per-state kernel discrepancy at most `epsilon_j`. Let the initial distributions differ by at most `epsilon_0`. A telescoping replacement argument and Markov-kernel contraction give

\[
\operatorname{TV}(P_0K_1\cdots K_r,\widetilde P_0\widetilde K_1\cdots\widetilde K_r)
\le \epsilon_0+\sum_{j=1}^r\epsilon_j.
\]

If the downstream ideal kernels have Dobrushin contraction coefficients `alpha_j <= 1`, use the telescoping order whose downstream factors are ideal to obtain

\[
\epsilon_0\prod_{j=1}^r\alpha_j+
\sum_{j=1}^r\epsilon_j\prod_{\ell=j+1}^r\alpha_\ell.
\]

The error bounds must be uniform over reachable interface states and interventions, or a separately justified state-restricted argument must replace this theorem. An isolated-system accuracy number is insufficient for an interacting composition.

## Preserving complementary value

An event probability, including classification error, changes by at most the target transcript TV discrepancy. If each of four risks defining `J = R_A + R_B - R_0 - R_AB` changes by at most `epsilon_target`, then

\[
J_{\rm actual}\ge J_{\rm ideal}-4\epsilon_{\rm target}.
\]

In the two-law construction, acquisition uncertainty alone yields the sharper `J >= 1/16 - delta_acquisition/4`. Adding target mismatch gives the sufficient positivity condition

\[
1/16-\delta_{\rm acquisition}/4-4\epsilon_{\rm target}>0.
\]

This is sufficient, not necessary. For an `r`-law chain, early marginal improvements are exponentially small in `r`; their robust certification requires correspondingly smaller errors. Ideal doubling does not remove this cost.

## Failure conditions

Independent-shot concentration fails for arbitrary correlations. A source that samples one bad calibration bit and repeats it forever has no improvement from repeated shots. Model membership also remains an assumption: context-dependent interactions can vanish on every calibration input and appear on unseen interventions. These limitations motivate the explicit nonlinear misspecification tests in the experiment suite.
