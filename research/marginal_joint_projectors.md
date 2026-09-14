# Joint half-filled and charged projector constraints

This work continues the exact six-site Hubbard certificate with a joint
weighted overlap bound. The new exactly accepted million-site lower density is
-1006233292669113/1562500000000000, approximately -0.6439893073082323.
The historical physical upper is -0.6106763470511881. Their interval width
is 0.03331296025704422 per site, 6.1323% narrower than the preceding interval.
The upper construction was referenced from its accepted receipt, not rerun here.

## Why separate penalties cannot help

The half-filled projector P_h and the sum P_c of four charged projectors
have disjoint spin-number supports. Therefore P_h+P_c <= I exactly.
Given K+k P_h+l P_c >= ell I with k,l >= 0, set delta=min(k,l).
The matrix for the reduced coefficients and lower value ell-delta is the
original matrix plus delta(I-P_h-P_c), so remains positive semidefinite.
With separate translated ceilings theta_h and theta_c, this changes the
certified density by delta(theta_h+theta_c-1)/5.
Here theta_h=0.5242147062 and theta_c=0.529758678, whose sum exceeds one.
Thus the independent two-penalty family is dominated by its boundary.
The charged-only boundary remains capped by the cyclic half-filled witness.

## Coupled construction and a dominance threshold

Let Q_r=P_h+r P_c, with r>0. A finite-window exact Gram bound
sum_{j=0}^{w-1} Q_{r,j} <= B I yields theta_J=B/w after averaging every
translation on a periodic chain longer than its support. The half-filled
bound may have a different window count: average that inequality separately.
There are no omitted projector windows in this periodic argument.

If alpha,beta >= 0 and

    K+(alpha+beta)P_h+beta*r*P_c >= ell I,

then the open-chain lower density is

    (ell-alpha*theta_h-beta*theta_J)/5 - (2*t+abs(V))/N.

The final correction removes the closing Hamiltonian bond after both
periodic inequalities have been combined.

For every r>0, local orthogonality gives

    Q_r <= r I+(1-r)P_h.

Write k'=alpha+(1-r)beta. Whenever k'>=0, the joint local inequality
implies the half-only inequality K+k'P_h >= (ell-r*beta) I.
The resulting half-only density minus the joint density is exactly

    beta * [theta_J-r-(1-r)*theta_h] / 5.

Consequently, **if** theta_J >= r+(1-r)*theta_h, that joint certificate
is dominated. For 0<r<=1 this applies to every alpha,beta>=0. This is a
conditional threshold on the supplied ceiling, not a claimed lower bound
on the optimal joint ceiling. A joint ceiling below the threshold is
necessary, but not sufficient, for a new energy improvement.

For r>1 and k'<0, the largest local penalty eigenvalue is beta*r.
Then K >= (ell-beta*r) I, and the unpenalized density exceeds the joint
density by [alpha*theta_h+beta*(theta_J-r)]/5 whenever that expression
is nonnegative. One must not multiply an upper bound on P_h by the
negative coefficient 1-r without reversing its direction.

## Exact weighted Gram replay

Embedded physical source columns C_j have integer amplitudes and squared
norms n_j. Put D_jj=1/n_j for the half source and r/n_j for each charged
source. For positive weights,

    C D C^T <= B I  iff  B D^(-1)-C^T C >= 0.

The accepting computation uses rational diagonal entries B*n_j/weight_j
and exact integer PSD checks. Numerical eigenvalue proposals must use the
symmetric matrix D^(1/2) C^T C D^(1/2). Feeding D C^T C to a symmetric
eigensolver is invalid; a draft discovery implementation doing this was
caught before any result was accepted.

Every exterior occupation and every spin block must be present. Odd source
embedding signs are constant on each exterior column, so they conjugate
the Gram by a diagonal sign matrix and do not change its spectrum.

## Verification status

The charged/joint backend passed its five focused standard-library tests.
The complete marginal regression then passed 474 tests in 331.122 seconds;
see results/marginal_graded_hubbard8/charged_projector_full_validation.log.
The combined energy path subsequently passed 476 regression tests in
463.443 seconds. The family-limit verifier added three focused tests, all
passing in 5.646 seconds; its final full regression passed479 tests in406.727 seconds.
The earlier malformed-source test used determinant63, which actually has
spin numbers(3,3); replacing it with the intended spin-(6,0) source corrected
that test expectation without changing the verifier.

The general physical-marginal cone, controlled accuracy-versus-cost scaling,
and transfer to generic molecular interactions remain unresolved.

## Accepted energy improvement

The new v4 energy path retains one rank-one update in each conserved local
spin sector. It replays the half bound and the joint bound, then proves
local positivity in all 94 reflection/spin blocks covering all 4,096 Fock
states. The largest local block remains200; no matrix cap was increased.

The tight r=1/2 joint ceiling B=595736407/250000000 is accepted within the
combined energy replay, giving theta_J=0.595736407. The separate overlap
receipt contains coarser accepted r=1/2,1,2 ceilings; do not confuse those
with the tighter numerical proposals. Independent binomial counts verify
all32 joint spin blocks, 1,280 total columns and maximum dimension132.

For reflected profiles with notation from the six-site report, the final
signed-density coefficients are

```
a=.147275    b=3.632038
p=.588477    q=1.181137
d=-.147865   e=.847147
alpha=.181963    beta=.094813    r=1/2
ell=-3.0680628
```

These give the exact periodic lower density -0.6439868073082323 and the
open-chain lower above after subtracting2.5/N. Signed density profile
coefficients were already admitted by the exact local verifier; onsite and
hopping profiles retain their required means and nonnegativity.

Freshly optimized nonnegative-profile half-only and joint certificates also
both accepted. Their periodic densities are -0.6442433633553013 and
-0.6440838534813463, respectively, so the joint penalty itself contributes
a verified improvement of about0.00015951 in that matched comparison.
Further signed-density tuning improves the joint result by about0.00009705.
The larger improvement over the old checkpoint includes ordinary profile
retuning and must not all be attributed to the new joint constraint.

Authoritative receipts are under
`results/marginal_graded_hubbard8/joint_projector/`:

* `profile_half_only_certificate_replay.json`
* `profile_joint_r1_2_certificate_replay.json`
* `signed_density/profile_joint_r1_2_certificate_replay.json`
* `independent_count_replay.json`
* `signed_density/combined_summary.json`

The combined standard-library replay took542.654 seconds in its observed
local run; this is not a controlled performance comparison. The numerical
search used a bounded SLSQP epigraph formulation with analytic eigenvector
gradients, correct affine derivative scales and fresh exact-CAR matrix
reconstruction before rounding. The earlier agent full-profile probe was
rejected for missing the derivative step divisor and falsely reporting a
fresh matrix check; its output is retained as a rejected diagnostic.

## The current relaxation's numerical limit is certified

Nine integer physical local vectors with positive rational mixture weights
form a trace-one local PSD matrix rho. Independent standard-library replay
reconstructs its physical CAR expectations. All six profile derivatives
vanish exactly, and both fidelities saturate exactly:

```
tr(rho P_h) = .5242147062
tr(rho (P_h + .5 P_c)) = .595736407
```

Taking rho's expectation in any local certificate and subtracting the
nonnegative penalty costs proves an upper limit on that certificate's
periodic lower density. This argument covers every reflected mean-correct
profile, including signed coefficients; it is not limited to the numerical
search box or a particular optimizer trajectory.

The exact rational cap is approximately -0.643986740577249. Together with
the accepting lower certificate it brackets the best attainable periodic
bound in this fixed-source, fixed-ratio, fixed-ceiling family within
6.67309833381364e-8 per site. See `signed_density/family_limit_certificate.json`
and `signed_density/family_limit_replay.json`; the exact cap fraction is in
those files. This is an upper limit on a relaxation's lower certificates,
not a physical variational upper bound on ground energy.

The mixture is compact (nine sources, roughly104KB JSON). Its exact
physical-expectation replay took2.224 seconds. It gives no fixed-cost
accuracy guarantee when the support, projector family or required accuracy
changes.

## A next physical constraint, actually violated

The limit mixture does not have matching left and right five-site marginals.
An exact diagonal entry of their difference at determinant102 is about
-0.008773379064595197. Its reflection partner is determinant612. Thus the
reflection-odd five-site diagonal operator

    Y = |102><102| - |612><612|

has exact expectation about -0.017546758129190394 in the marginal difference.
The six-site correction Y_left-Y_right is reflection-even and telescopes to
zero on the periodic chain. This identifies physical information absent
from the just-capped relaxation. The exact rational separation is recorded
in `signed_density/overlap_consistency_separator.json`.

A first numerical optimization with only this one extra correction returned
zero correction and no improvement. Violating one chosen dual mixture does
not prove that every dual optimum violates that same constraint. Jointly
imposing several overlap constraints, or changing the weighted ratio/source
family, remains necessary to test. The exploratory telescoping proposal is
not an accepted energy-path extension.
