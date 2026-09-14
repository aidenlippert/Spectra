# Review of estimated signed-moment planning bounds

Assume each hidden model `m` has signed label `s_m` and prior `p_m`. For an
action `a`, let `L_{ma}=Pr(+|m,a)`. A two-shot policy has four outcome leaves;
its signed success margin can be expanded into the raw moments

`mu0 = sum_m p_m s_m`, `mu_a = sum_m p_m s_m L_{ma}`, and
`G_ab = sum_m p_m s_m L_{ma}L_{mb}`.

For a fixed policy, if every estimate has uniform errors
`|dmu0|<=e0`, `|dmu_a|<=e1`, and `|dG_ab|<=e2`, the leaf expansion gives

`|d margin| <= e0 + 4 e1 + 4 e2`.

Therefore the fixed-policy risk error is at most

`kappa = e0/2 + 2e1 + 2e2`.

If the planner maximizes estimated margin over the finite policy class, the
selected policy's true regret against the best policy is at most `2kappa`:

`regret <= e0 + 4e1 + 4e2`.

The factor two is necessary: one error bound applies to the selected policy
and one to the true optimum. This is a finite robust-planning statement, not a
claim that estimated moments recover the physical model.

## Sampling bound

With labelled calibration episodes sampled from the same hidden-model prior,
each signed contribution lies in `[-1,1]`. Hoeffding plus a union bound gives

`N >= (1/(2 e^2)) log(2K/delta)`

for simultaneous additive error `e` over

`K = 1 + A + A(A+1)/2`

moments, counting the symmetric Gram entries. If moments are estimated with
different tolerances, use the corresponding union bound rather than silently
substituting one common `e`. Calibration labels, action costs, and the number
of episodes must be charged separately from policy execution.

## Required caveats

The empirical moment tuple need not be physically realizable. It may violate
positivity, symmetry, or probability constraints even though every component
estimate is individually plausible. The regret theorem still applies to the
finite linear surrogate if the same estimates are used consistently, but a
physical-model certificate must additionally project onto, or produce an
interval certificate for, the realizable moment polytope. Projection error must
be added to `e0,e1,e2`.

Zero-probability branches create no contribution to risk and need no leaf
decision. A planner must represent them explicitly and avoid dividing by an
estimated branch probability. Repeated actions are included by allowing
`a=b`; then `G_aa` is required and counted once in the symmetric table. If
episodes are adaptively selected, the iid Hoeffding statement no longer applies
without a martingale argument and a declared sampling policy.

The bound assumes conditional independence of the two readouts given the same
hidden model. Correlated detector noise, drift between shots, or a changed
latent state invalidates the product `L_{ma}L_{mb}`; their contribution must be
represented by additional moments and charged in the error budget.

This result supports a bridge from labelled source data to certified bounded
planning. It does not establish cumulative scientific capability: positive
transfer to a fresh task and discovery acquisition cost remain separate tests.
