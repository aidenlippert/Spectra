# Next bounded construction test: inexpensive overlap of ordinary partitions

The V11 normalized-Frobenius rejection guard is now the strengthened conventional
Taylor baseline. The exact matrix route lost and is closed. No acquired method
exists, and no reserved evaluation has been opened.

## Concrete existing calculation

V8 residual row10 is n3 XXZ, gamma2, T1/2, degree11. Its weighted ordinary
partition gives error0.00101656394198456, just above the required0.001.
The next Taylor degree passes. The old full overlap proposer used69 candidate
groups,16 sweeps and large independently rationalized coefficients; it reduced
the bound but increased complete cost. That implementation remains rejected.

The agent's initial conclusion that no small family or one-pass split was
supported was too strong. It also demanded matching the old0.000965 bound,
although the requirement is only0.001. Root tested actual certificates:

| Group family | Sweeps | Integrated bound | Passes .001? | Norm construction + checks |
| --- | ---: | ---: | --- | ---: |
| Full original family | 1 | .000989288 | yes | .01668 s |
| Two ordinary partitions | 0 | .001016564 | no | .00313 s |
| Two ordinary partitions | 1 | .000998348 | yes | .00372 s |
| Two ordinary partitions | 4 | .000986802 | yes | .00382 s |
| Two ordinary partitions | 16 | .000986012 | yes | .00418 s |

The one-sweep two-partition witness has14 groups, not61. The union is obtained
from the existing weighted and firstfit partitions. These were norm-only
diagnostics on an already public development residual, not end-to-end timing
headroom or method acquisition. Source hash, costs and checker receipts are in
`results/v11/cover_schedule_probe.json`. The uncorrected agent conclusion is
preserved under `research/v11/invalidated/agent_cover_structure.md`.

## A distinct remaining cost mechanism

The old proposer converts each optimizer coordinate independently with
`Fraction.from_float(...).limit_denominator(10**9)`, then exactly reconciles
each coefficient. Incommensurate denominators create much larger denominators
in the reconciliation sum. This is avoidable without weakening a certificate.

After one cyclic coordinate sweep on only the two ordinary partitions, express
each incidence as a nonnegative weight of its original coefficient. Quantize
weights to a common dyadic denominator M=2^20; allocate floors to all but one
incidence and the exact remaining integer weight to the last. Enforce nonnegative
integer weights summing exactly to M; untrusted floating proposals can be clipped
or refused. Then each split is `c_P a_GP/M`, and exact reconstruction follows
from `sum_G a_GP=M`. Use the unchanged V8 overlap checker for anticommutation,
root bounds and exact reconstruction. It need not trust convergence or the
optimizer's norm estimate. This is conventional rational representation, not
novel physics or an acquired method.

No speedup is assumed. Only one degree can be saved on this case, so the saved
construction/checking work may be smaller than the extra cover cost. The first
test must measure the entire adaptive calculation against V11 Frobenius-guarded
Taylor, reusing already computed ordinary partitions only with their cost still
counted. Retain the old physical input, tolerance and support caps. Include the
cover's construction, failed proposals, and a fresh complete evolution check.
If it loses, close it; do not rescue it by comparing with the now-weaker V8
baseline or reporting only the smaller residual norm.

A pass would only permit development expansion and proper paired timings.
It would not license a compounding claim. Any subsequent method learner must
beat this strengthened conventional library and still pass the m1-to-m2 causal
acquisition and two-heldout-increment requirements.
