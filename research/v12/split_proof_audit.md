# V12 split proof audit

This audit covers the proposed one-coordinate-sweep cover made from the union
of the existing weighted and firstfit ordinary partitions. It uses the stated
14-group witness and the common dyadic denominator

\[
 M=2^{20}=1,048,576.
\]

The reported one-sweep integrated norm is `0.000998348`, which is below
`0.001` as a norm-only result. It does not establish the required complete
calculation or timing headroom. The construction remains a conventional exact
representation of a split certificate.

## Preconditions and scope

For each nonzero real input coefficient `c_P`, the optimizer supplies
nonnegative candidate weights `w_{G,P}` over the incidences of `P`. The
candidate is untrusted. The trusted construction must use only finite,
nonnegative values after validation/clipping and must normalize them so that

\[
 \sum_{G\ni P} w_{G,P}=1.
\]

If an input coefficient is zero, non-real, non-finite, or missing from the
incidence union, the proposal is refused. The signed coefficient is then
split as `c_P a_{G,P}/M`; the sign is inherited from `c_P`, and all arithmetic
used for the final witness is exact integer/rational arithmetic. There is no
trusted use of a floating objective, convergence claim, or optimizer output.

The union is deduplicated by group label. Every group must pass the unchanged
V8 anticommutation check, and every Pauli must have at most two retained
incidences unless the construction explicitly rejects the witness. The
ordinary partition groups are allowed to have one incidence; a one-incidence
coefficient simply receives the full weight `M`.

## Integer allocation invariant

Let `k_P` be the number of retained incidences of `P`, with `1 <= k_P <= 2` for
this audit. Select one incidence as the largest-incidence destination. For all
other incidences set

\[
 a_{G,P}=\lfloor M w_{G,P}\rfloor,
\]

and set the destination to the exact remaining integer

\[
 a_{G_*,P}=M-\sum_{G\ne G_*}a_{G,P}.
\]

Because each floor is between zero and its argument and the normalized weights
sum to one, the remainder is an integer at least zero and at most `M`. Thus all
integer weights are nonnegative and sum exactly to `M`; no floating remainder
is admitted. If the implementation requires every retained incidence to be
strictly positive, it must either seed each incidence with one unit before
flooring or drop zero-weight incidences and renormalize. The simpler
floor-plus-remainder rule is already sufficient for the V8 checker because
zero split coefficients are omitted, but it must never emit a zero incidence
to `check_fractional_cover` (that checker rejects zero coefficients).

For a nonzero `c_P`, `c_P a_{G,P}/M` is nonzero exactly when `a_{G,P}>0`, and

\[
 \sum_{G\ni P} c_Pa_{G,P}/M=c_P.
\]

Therefore reconstruction is exact independently for every Pauli. Summing over
Paulis gives exact operator reconstruction. The residual branch must use the
coefficient `-G c_m` (or its exact rational equivalent): negating every split
coefficient preserves both reconstruction and each group square norm. Using
`+G c_m` would certify the wrong residual despite identical norm values.

## Quantization inflation bound

Write `q_{G,P}=a_{G,P}/M` and `e_{G,P}=q_{G,P}-w_{G,P}`. For every floored
non-destination incidence, `-1/M < e_{G,P} <= 0`; the destination error is the
negative sum of those errors. Consequently

\[
 \sum_{G\ni P}|e_{G,P}|
 \le {2(k_P-1)\over M}.
\]

Let `v_G` be the unquantized coefficient vector in group `G` and `d_G` the
quantization perturbation. The triangle inequality gives the rigorous bound

\[
 \sum_G\|v_G+d_G\|_2
 \le \sum_G\|v_G\|_2
   + {2\over M}\sum_P |c_P|(k_P-1).
\]

With `k_P<=2`, the additive inflation is at most

\[
 {2\over 2^{20}}\sum_P|c_P|
 = {\|c\|_1\over524288}.
\]

This is an upper bound, not an observed error estimate. If the construction
uses outward exact square-root intervals as in `_witness`, their endpoint
rounding is an additional checker bound and must be included in the claimed
bound; the float objective cannot substitute for it. For a Taylor residual
weighted by `factor`, the integrated inflation is bounded by `factor` times
the displayed additive term. The available margin from `0.000998348` to
`0.001` is `0.000001652`; it is therefore necessary to compute the exact
residual-specific bound and leave room for root-interval outward rounding and
all other construction costs. The generic `||c||_1/524288` bound alone does
not prove that the margin is sufficient.

## Checker and cost implications

The unchanged `check_fractional_cover` already checks group anticommutation,
nonzero coefficients, exact rational reconstruction, each root upper bound,
and the claimed sum. `check_evolution_cover` then recomputes residual records,
checks the horizon and witness labels, and compares the exact weighted total
to tolerance. These remain the trusted rules; no new checker rule or trusted
optimizer state is needed.

The 14 groups reduce group count relative to the rejected large overlap
proposal, but each incidence still costs an exact coefficient addition and
square, and each within-group pair costs an anticommutation check. Dyadic
weights require at most 20 denominator bits plus the numerator range needed
for `M|c_P|` in exact rational arithmetic. The denominator does not multiply
across incidences, unlike independently limited denominators. Nevertheless,
the complete cost must include partition construction/reuse, failed proposals,
integer allocation, exact witness checking, and a fresh complete evolution
check. No speed or headroom conclusion follows from the norm-only receipt.

## Audit disposition

The split proof is sound under the stated validation gates: finite real signed
coefficients, nonnegative normalized weights, exact floor/remainder allocation,
deduplicated union, and unchanged V8 checks. The proof establishes exact
reconstruction and a conservative dyadic norm-inflation bound. It does not
establish a passing V12 result, because the residual-specific inflation,
outward root intervals, and complete acquisition-inclusive timing have not
been supplied. The one-sweep `0.000998348` result should therefore remain a
development diagnostic until a fresh `check_evolution_cover` receipt and full
cost comparison are available.
