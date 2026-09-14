# Trustworthy prediction under adaptive physical design

## Claim that is actually supportable

Let `D={1,...,M}` be a finite set of physically executable design/policy
choices. Design `d` has bounded loss `J(d) in [0,1]`. At round `t`, an
adaptive policy chooses `d_t` from all previous observations and receives
`Y_t = J(d_t)+xi_t`, where `E[xi_t | F_{t-1},d_t]=0` and `xi_t` is conditionally
`sigma`-sub-Gaussian. Suppose the procedure maintains intervals
`[L_t(d),U_t(d)]` satisfying the simultaneous event

```text
P( for every t <= B and every d in D: J(d) in [L_t(d),U_t(d)] ) >= 1-delta.
```

This is the needed assumption: it is uniform over the adaptive stopping time,
the adaptively selected designs, and the finite budget. It can be obtained by
time-uniform confidence sequences (Howard et al., 2021), with a union bound
over `M` designs, or by a separately justified adaptive-data-analysis
mechanism. It is not implied by an ordinary random train/test error.

If the returned design `d_hat` satisfies
`U_B(d_hat) <= min_d U_B(d) + eta` for minimization, then on the simultaneous
event

```text
J(d_hat) <= min_d J(d) + 2 r_B + eta,
```

whenever every interval has radius at most `r_B`. For constraints `g_j(d)<=0`,
accepting only when `U^g_{j,B}(d)<=0` gives physical feasibility on the same
event. A constrained optimizer must compare only against a competitor already
certified feasible (or report that no such competitor is certified); minimizing
an upper bound over all designs can select an uncertified infeasible competitor.
Thus a finite-budget adaptive design guarantee is possible for a
declared finite design class, an explicit stochastic observation model, and a
uniform confidence mechanism.

### Proof

Let `d*` minimize `J`. On the simultaneous event,
`J(d_hat) <= U_B(d_hat) <= U_B(d*)+eta <= J(d*)+2r_B+eta`.
For a constraint, `g_j(d)<=U^g_{j,B}(d)<=0`. The argument is purely deterministic
after the uniform event is established. If the minimizer does not exist, use an
arbitrarily near-minimizing sequence.

Under the stated mgf convention `E exp(lambda xi) <= exp(lambda^2 sigma^2/2)`,
a fixed-`n` sub-Gaussian interval has radius
`sigma sqrt(2 log(2MB/delta)/n)`. For countably many sites `i`, functions `j`,
and sample sizes `n>=1`, an explicit union allocation is

```text
alpha_{i,j,n} = 6 alpha / (pi^2 N K n^2),
r_j(n) = sigma_j sqrt((2/n) log(pi^2 N K n^2/(3 alpha))).
```

The corresponding two-sided bounds hold simultaneously with probability at
least `1-alpha` (using `sum_n 6/(pi^2 n^2)=1`). If allocation is predictable
and observations have martingale sub-Gaussian increments, the same statement
follows from the exponential supermartingale and optional-stopping validity;
confidence-sequence boundaries are the sharper general replacement. The exact
width still depends on the sampling design and variance method.

## What this does and does not certify

There are three distinct errors. Numerical error is error in evaluating a
declared model; it can be bounded by a discretization, truncation, or residual
certificate. Statistical noise is variation in observations conditional on the
design; confidence sequences address it under their martingale/noise
assumptions. Model-class error is the gap between the declared `J` and the
physical response; no statistical interval for a misspecified `J` bounds that
gap. The theorem certifies the physical response only if `J(d)` denotes that
response and the observation model is valid. A surrogate interval around a
simulator certifies simulator accuracy, not physical applicability.

The finite class is a deliberate boundary. For a continuous or history-rich
space, replace the union bound by a uniform complexity condition (for example,
a proven RKHS norm/kernel model, Lipschitz covering bound, or a reusable
holdout mechanism). Those assumptions are substantive physical claims and must
be declared and checked; GP or Bayesian posterior width alone is not a
frequentist physical guarantee unless its coverage assumptions are proved.

## Elementary obstruction: average accuracy can fail exactly at the design

For each `M`, let the design space be `{1,...,M}` and let the true losses be
`J(d)=0` for every `d`. Define a surrogate `J_hat(1)=-1` and
`J_hat(d)=0` otherwise. Under the uniform design distribution its absolute
error is `1/M`, which tends to zero, yet optimization always selects `d=1`,
where the claimed improvement is fictitious. A variant with `J(1)=1` and
`J_hat(1)=-1` gives a two-unit failure while average absolute error is `2/M`.

This is not a pathology of optimization implementation; it is an information
obstruction. Any guarantee based only on average held-out error permits an
unseen or low-probability region containing the optimizer's selected design.
Adaptive querying makes the selection event part of the guarantee. Dwork et al.
show the corresponding statistical issue for adaptively chosen analyses: naive
reuse of a fixed holdout can lose validity, while reusable-holdout mechanisms
preserve it under explicit stability/privacy conditions.

## Primary sources and relevance

* Howard, Ramdas, McAuliffe, and Sekhon, “Time-uniform, nonparametric,
  nonasymptotic confidence sequences,” *Annals of Statistics* 49 (2021),
  [author paper](https://arxiv.org/abs/1810.08240).
  Confidence sequences are uniform over time and support optional stopping;
  they supply the statistical part of the simultaneous event.
* Dwork, Feldman, Hardt, Pitassi, Reingold, and Roth, “Preserving Statistical
  Validity in Adaptive Data Analysis,” STOC (2015),
  [paper](https://www.cs.toronto.edu/~toni/Papers/adaptive-stoc.pdf).
  It formalizes why adaptively selected queries require mechanisms beyond
  ordinary fixed-sample generalization.
* Sui, Gotovos, Burdick, and Krause, “Safe exploration for optimization with
  Gaussian processes,” ICML 2015, [PMLR paper](https://proceedings.mlr.press/v37/sui15.html).
  It is a primary example of adaptive optimization with high-probability safety
  under explicit GP regularity and confidence assumptions; those assumptions
  delimit applicability rather than establish universal physical validity.

The honest mission-level target is therefore: construct experiments and models
that maintain a uniform, physically grounded confidence event over the adaptive
designs they will actually propose, and explicitly detect/refuse when numerical
error, statistical noise, or model-class applicability is not bounded. Without
that event, “trustworthy prediction under invention” is not established.
