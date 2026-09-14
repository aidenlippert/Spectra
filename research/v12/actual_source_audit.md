# V12 actual source audit

This is a read-only audit of `experiments/v12_small_cover.py`,
`experiments/v12_small_cover_taylor.py`, and `tests/test_v12_small_cover.py`.
The implementation and tests support exact certificate soundness. They do not
make the floating proposal trustworthy, and they do not by themselves establish
a general normalized-weight theorem or a cost/headroom win.

## Exact reconstruction

`small_cover` first applies `clean(op, n, 40)`. Thus zero input coefficients are
removed rather than rejected. An all-zero input returns an empty witness. For a
nonempty cleaned operator, each label is required to occur in exactly one group
of each supplied ordinary partition; the union is deduplicated by sorted group
labels. Every surviving label is required to have one or two incidences.

For each incidence, the code computes finite floating scores from the sweep,
takes absolute values, and normalizes by their floating sum. If that sum is
zero, it uses equal fallback scores. It then exports integer allocations with
`M=1<<grid_bits`: each non-anchor allocation is a clamped integer floor of its
score fraction, and the anchor receives the remaining integer. `remaining`
starts at `M` and only decreases by a nonnegative allocation no larger than it,
so every allocation is an integer in `[0,M]` and allocations for each Pauli sum
exactly to `M`. The emitted coefficient is `op[p] * F(a,M)`, preserving the
input sign, and `_witness` omits zero entries.

Therefore the exact reconstruction property is independent of the quality,
normalization, or convergence of the float sweep: for every retained nonzero
`p`, the emitted coefficients sum exactly to `op[p]`. The unchanged
`check_fractional_cover` additionally checks the exact rational sum, group
anticommutation, nonzero emitted incidences, root upper bounds, and claimed
bound. The tests exercise both reconstruction and rejection after changing an
upper bound or coefficient.

The checker does not rely on a normalized-weight theorem. It sees only the
exported exact coefficients. The implementation's score normalization is an
internal route to choosing integers, not a trusted premise.

## Float influence and the optional estimate

Floats influence only the proposed allocation: initial equal shares, the
incremental square/roots in the coordinate sweep, score ordering, and the
integer floors. Nonfinite values are rejected. The exact witness export and
the unchanged checker remove float rounding and cancellation from the trusted
certificate. A bad finite proposal can still produce a looser exact bound; it
cannot silently produce a false accepted reconstruction.

The earlier idealized description said “floor all but largest normalized
weights.” The source implements the equivalent safety shape using absolute
scores and sequential floors, with the anchor selected by maximum score. It
does not prove that the anchor is the mathematically largest ideal weight, and
no such fact is required for exactness. Likewise, the generic inflation bound
`2/M * sum |c_P|(k_P-1)` is only an optional estimate when `w` denotes a
specified exact normalized target and each exported integer is its floor-plus-
remainder quantization. It is not a source-level correctness condition and
must not be presented as a theorem about the arbitrary float scores used here.
The required guard is the exact checker result and its exact summed root bound.

## Taylor hook and guard completeness

`small_cover_taylor` keeps the V11 integer recurrence and denominator/bit
budgets. At each order it computes the existing `norm_rejection` guard before
building ordinary partitions. A rejected order is recorded and skipped. For a
non-rejected order, it computes all three ordinary probes (`l1`, `firstfit`,
`weighted`) and returns immediately on an ordinary success. The small-cover hook
is reached only when the derivative has at most 40 terms and the best ordinary
bound is within `11/10` of tolerance. It uses exactly the weighted and firstfit
partitions, performs one sweep, and records `cover_calls` including construction
work.

On a cover success, `evolution_witness` inserts the cover only at the requested
residual degree and negates every exported coefficient. This is the required
residual sign `-G c_m`; the test deliberately negates it again and confirms
that `check_evolution_cover` rejects the resulting mismatch. The evolution
checker independently recomputes residual records, checks labels and horizon,
rechecks each cover exactly, and compares the exact weighted total to
tolerance. The complete lower-degree test passes this route at degree 11 and
checks that the ordinary guarded baseline needs degree 12.

The hook's `weighted_bound` is an untrusted caller-supplied claim at construction
time, but it is not accepted on faith: `check_evolution_cover` recomputes the
exact total from the witnesses and requires equality with the claimed bound.
The implementation therefore has the needed guard/check path for correctness.

## Remaining limits

The tests establish exact reconstruction, rejection behavior, residual sign,
disabled-hook equivalence, and the fixed development case. They do not establish
that the complete adaptive calculation wins on acquisition-inclusive cost or
has sufficient timing headroom. The 14-group/one-sweep norm result is a
development result until fresh full receipts include construction, failed
proposals, ordinary partition work, exact checks, and the complete evolution
check. No novelty, convergence, or matching of any older bound follows from
this source audit.
