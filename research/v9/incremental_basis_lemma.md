# Incremental projected-basis extension lemma

This is a bounded conventional-computation proposal for the next projected
Taylor construction.  It does not constitute a learned method or a held-out
result.

Let `B ⊂ B'` be coordinate subspaces with projectors `P_B` and `P_B'`, and
let `G` be the exact sparse generator action.  Embed the old coefficients in
`B'` and define

\[
 c_0'=c_0,\qquad c_{k+1}'={P_{B'}G c_k'\over k+1},
\]

while the old projected recurrence is

\[
 c_0=c_0,\qquad c_{k+1}={P_BG c_k\over k+1}.
\]

The inclusion makes `P_B'` act as the identity on every vector already in
`B`.  Assume that for `k<r`, the newly available component vanishes:

\[
 (P_{B'}-P_B)G c_k=0.
\]

Then `c_k'=c_k` for every `k≤r`.  The proof is induction.  The base case is
the shared initial coefficient.  If `c_k'=c_k` and the displayed leakage
condition holds, then

\[
 c_{k+1}'=P_{B'}Gc_k/(k+1)=P_BGc_k/(k+1)=c_{k+1}.
\]

The first coefficient that can differ is therefore the coefficient after the
first nonzero leakage order.  No numerical tolerance is involved: the test is
exact support membership and exact rational zero.

## Forced-difference recurrence

Set `δ_k=c_k'-c_k` and `Q=P_{B'}-P_B`.  Subtracting the two recurrences gives

\[
 \boxed{\quad
 δ_0=0,\qquad
 δ_{k+1}={P_{B'}Gδ_k+QGc_k\over k+1}.\quad}
\]

This identity separates propagation of an already introduced difference from
the new forcing supplied by the old trajectory.  It also shows exactly what
must be computed at an extension: the old raw derivative `Gc_k` is sufficient
for the forcing term, while the new projected action is needed only on `δ_k`
and on the newly generated suffix.  If the forcing term is zero through
order `r-1`, the induction above gives `δ_0,...,δ_r=0`; at order `r`, the
first nonzero `QGc_r` determines the change.

## What can be reused

For a fixed polynomial order `m`, a restart from `B` to `B'` can reuse the
following exact proposer-side work through the unchanged prefix:

* the coefficients `c_0,...,c_r`;
* cached raw derivatives `Gc_0,...,Gc_{r-1}` and their sparse rational maps;
* residual coefficients whose formulas use only those unchanged coefficients
  and derivatives; and
* norm witnesses for residual operators that are literally identical.

The suffix must be generated from the first affected order onward.  In
particular, a cached `Gc_k` is reusable only when its input coefficient is
unchanged; the old derivative cache cannot be applied to a changed `c_k'` by
label alone.  Extending polynomial order from `m` to `m+q` is the same
recurrence with a longer suffix: all coefficients and raw derivatives through
`m` remain valid, but each new coefficient and residual must be generated and
witnessed.

For a continuous single piece, an unchanged prefix of residual records can be
copied into a proposer transcript because the operator and integration weight
are identical.  For piecewise trajectories, this reuse is valid only within a
piece whose initial value, duration, and preceding interface value are
unchanged; a changed endpoint can alter the next piece's jump and all of its
coefficients.

These are construction savings, not permission for the independent checker to
trust a cached certificate.  A fresh checker must recompute the candidate's
residuals from the submitted coefficients, verify every reused norm witness
against those recomputed operators, and check the complete claimed integral,
jumps, horizon, and tolerance.  It may use a separately implemented exact
prefix hash or equality check as an optimization only after the same data are
available to the checker; an unverified proposer cache is not evidence.

## Cost condition

Let `C_k` denote the charged cost of applying `G` to a coefficient at order
`k`, including sparse rational multiply-adds, and let `W_k` denote its norm
witness cost.  A restart with unchanged prefix through `r` saves roughly

\[
 \sum_{k<r}(C_k+W_k)
\]

of proposer work, while adding the cost of exact leakage tests and the suffix
from `r` onward.  If the basis diagnostic is run every round and `r` is small
(the new basis direction is immediately active), this saves little and can be
slower because each round pays discovery, basis insertion, duplicate checks,
and cache bookkeeping.  If several rounds produce a long unchanged prefix,
the saving is material, especially when raw generator applications dominate.
The independent checker still pays its normal full replay unless a separately
audited incremental checker is introduced, so proposer savings alone cannot
establish complete-cost headroom.

The meaningful preflight is therefore a development comparison over the
existing `n=3,4,6`, degree-`≤24` grid: report prefix length, diagnostic cost,
raw-derivative reuse, suffix cost, witness cost, fresh-check cost, and total
wall time for every restart.  A useful case must beat the integer-Taylor or
current projected baseline after all restarts and failed candidates are
charged.  The result would establish an incremental implementation identity,
not a new acquired construction rule.
