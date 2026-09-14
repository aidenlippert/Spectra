# Review of the scalar-probe Gaussian lower bound

The proposed construction is a useful finite obstruction to one-shot scalar
readout. In orthonormal coordinates, let `c=3/5`, `s=4/5`,

`u_a=c e1+a s e2`, `t_a=-a s e1+c e2`,

and `v_ab=c t_a+b s e3`, with independent signs `a,b`. The vectors `u_a,t_a`
are orthonormal for each `a`. A scalar probe `x` has signal coefficients
`x·u_a` and `x·t_a`. Summing over `a` gives

`sum_a (x·u_a)^2 >= (9/16) sum_a (x·t_a)^2`,

Indeed, after multiplying by two this is
`2(c^2 x1^2+s^2 x2^2) >= (9/16)2(s^2 x1^2+c^2 x2^2)`.
Thus at least one nuisance mode satisfies
`(x·u_a)^2 >= (9/16)(x·t_a)^2`.

For that mode, the two covariance projections have midpoint `D` and
half-difference `E` obeying

`D >= Lambda*(9/16)(x·t_a)^2 + lambda*s^2(x·e3)^2`,
`|E| = 2 lambda*c*s |x·t_a||x·e3|`.

Using `2pq <= p^2+q^2` yields

`|E|/D <= (4/5) sqrt(lambda/Lambda) = r`.

The bound is scale invariant and remains valid when the nuisance-independent
PSD baseline `B` is added, since `B >= I` only increases the midpoint.

For centered one-dimensional Gaussians with variances `D+E` and `D-E`, an
explicit divergence bound is required. With `rho=|E|/D<1`, put
`t=2rho/(1-rho)`. The variance ratio is `1+t`; using
`log(1+t)>=t-t^2/2` gives `KL<=t^2/4`. Pinsker therefore gives
`TV<=t/(2sqrt(2))<=rho/(1-rho)`. The simpler `rho/(1-rho)` bound is safe but
must be stated as a deliberately loose inequality and checked for `rho<1`.
Consequently the nuisance mixture has TV at most
`0.5*(1+r/(1-r))`, and equal-prior target Bayes risk is at least

`1/4 - r/[4(1-r)]`.

For `Lambda/lambda >= 10000`, `r <= 1/125`, giving `risk >= 123/496`.
This lower bound covers arbitrary scalar probes and arbitrary postprocessing,
because postprocessing cannot increase total variation.

The warm construction chooses `q_a=-s t_a+c e3`. It nulls `u_a`, and its
projection on `v_{a,+}` is zero while its projection on `v_{a,-}` has magnitude
`2cs=24/25`. Thus it converts the target sign into a variance change. A fixed
threshold such as `|Y|>20` is only valid after an explicit scale choice and
tail-error calculation; `B>=I` alone does not make that threshold reliable.
The certificate should instead specify the likelihood, calibrated variance,
and target error, then use the exact Gaussian tail or a rational conservative
bound.

## Expected adaptive query lower bound

Let `r0=123/496-10^-9` be the one-query Bayes-error floor and let the target
prior be balanced. For a randomized variable-horizon learner targeting
`delta=1/20`, write `x=P(N>=1)` and `y=P(N>=2)`. Zero-query branches have risk
at least `1/2`, hence `x>=1-2delta`. Starting after one query, the risk is at
least `r0`; continuing can reduce posterior error by at most `1/2` on the
`N>=2` event. Therefore

`R >= .5(1-x)+r0*x-.5*y`,

so `E[N]>=x+y >= 1-2delta+2r0*x >= (1-2delta)(1+2r0)`, approximately `1.34637`.
The argument permits randomized stopping, arbitrary controls, and unlimited
potential query depth; it uses only the independent hidden mode and the
one-query TV floor. It also applies after conditioning on a retained nuisance
mode. With a fixed acquisition charge of one and four warm target queries,
the warm sequence costs five, versus at least four times the cold lower bound
(`about 5.38548`) across four fresh stages. The comparison is a net expected
query claim with every physical reset charged equally; a stateful conventional
learner with the same retained mode is a fair tie.

The extension to a third nuisance mode is valid by adding an orthogonal
coordinate plane and conditioning on previously identified modes. Its claim
is conditional: if the nuisance mode changes between stages, the warm probe
need not remain aligned and the transfer advantage disappears.
