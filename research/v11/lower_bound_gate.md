# V11 lower-bound feasibility gate

This is a conventional rejection-only optimization. It is not an m1 result,
novelty claim, or change to the V8 baseline. The gate may skip proposal work
only when the exact tail lower bound already exceeds the tolerance.

## Audit of the guarded implementation

The implementation in `experiments/v11_guarded_taylor.py` uses
`threshold=floor(right/left)` for the max guard and
`threshold=floor(right^2/left^2)` for the Frobenius guard. Since the inspected
values and the accumulated square sum are integers, `value > floor(x)` is
equivalent to `value > x` under the strict comparison. Thus equality proceeds,
and the comparator itself preserves the first accepted polynomial and grouping
choice. The derivative denominator is correctly
`q0*q^(k+1)*k!`; the Taylor coefficient fraction is still constructed before
the gate, while derivative fractions and all three norm proposals are skipped
only after a rejection. The final `max_order` iteration still computes the
same next integer derivative as V8 and then raises if no certificate exists.

The guarded caller establishes the relevant bounds before this internal helper:
each derivative numerator has at most 8192 bits, there are at most
`gen.max_terms <= 4096` entries, and the helper checks `left` and `right` at
most 16384 bits. Therefore `left^2` and `right^2` are bounded by 32768 bits;
each `value^2` is bounded by 16384 bits; and the running square sum is bounded
by 16396 bits (the extra 12 bits cover at most 4096 terms). The intended cap
semantics should document this distinction: 16384 bits for unsquared
comparison operands and 32768 bits for squared comparison temporaries (with
the corresponding 16396-bit sum bound). `max_comparison_bits` records the
unsquared operand width, not the largest temporary. Under the validated public
caller these are a priori bounds, so the absence of repeated checks inside
`norm_rejection` is not an actual counterexample or correctness bug. The helper
is internal, not an untrusted certificate boundary; malformed direct calls
would be outside its contract.

## Bound

At Taylor order `k`, let the exact derivative tail be

```text
R = sum_p r_p P,
r_p = b_(k+1,p) / D_k,
D_k = q0 * q^(k+1) * k!.
```

The V8 power-basis integration weight is
`w_k = T^(k+1)/(k+1)`. For any valid norm witness, each anticommuting group
has an upper bound at least its exact normalized-HS norm. Therefore

```text
sum_groups upper(group) >= sqrt(sum_p r_p^2) >= max_p |r_p|.
```

Consequently, if

```text
w_k^2 * sum_p r_p^2 > epsilon^2,                         (1)
```

no witness accepted by the V7 checker can pass at this order. The weaker
`w_k^2*max_p(r_p^2)>epsilon^2` test is also sound, but the sum-of-squares test
is the useful default.

## Exact integer comparison

Represent `w_k=Wn/Wd` and `epsilon=En/Ed` as reduced nonnegative integer
fractions, and compute `S=sum_p b_(k+1,p)^2` using unbounded integers. Test
(1) without constructing coefficient fractions:

```text
Wn^2 * S * Ed^2 > En^2 * Wd^2 * D_k^2.                    (2)
```

The comparison is strict. Equality does not reject: a valid witness could
touch the tolerance exactly, subject to the checker’s exact claimed-bound
comparison. With `T=a/t`, one may use unreduced positive powers in (2), but
reducing `Wn/Wd` first limits intermediate growth. Cancel common factors
between the two sides before multiplication when this is useful; never use
floating point or a rounded square root.

For `epsilon=0`, (2) rejects exactly when `S>0` and `w_k>0`; a zero tail is
accepted and proceeds through the ordinary empty-witness path. If `T=0` is
disallowed by the existing duration contract, retain that contract rather than
special-casing it here.

## Preservation properties

The gate can only reject an order whose best possible checker bound is already
over tolerance. If it does not reject, all three existing V8 choices
(`l1`, `firstfit`, and `weighted`) run unchanged. Hence it preserves the first
accepted order and the norm-choice result exactly: it removes no feasible
candidate and changes no grouping, witness, or bound. At the strict boundary,
the `>` test is required for this statement; replacing it with `>=` could skip
a valid exact certificate.

The gate must run before each expensive `norm_witness` call, after computing
the same integer derivative numerator `b_(k+1)` used by V8. It must not alter
the derivative, denominator, order cap, or probe ordering. A rejected order
should produce an explicit rejected/ infeasible probe record (with the exact
integer comparison result) so the accounting remains auditable; it must not be
reported as a successful certificate.

## Caps and cost accounting

The sum of squares and both sides of (2) are temporary exact integers. Apply
the existing rational/integer bit caps before allocating or multiplying them;
if a cap would be exceeded, use the established arbitrary-integer refusal or
budget failure path. Do not silently convert to a machine integer. Count the
`S` accumulation, integer squares, exact scaling/comparison operations, and
the order probe in proposal cost. Skipped sorting, grouping comparisons,
square enclosures, and witness construction are genuine savings, but the gate
cost and all already-computed derivative/action work remain charged.

The gate is part of the same immutable-checking budget semantics: it cannot
waive endpoint conversion, fresh checking, term/support caps, rational bit
limits, or order limits. It only establishes impossibility for that order.

## Methodological status

This is the standard Euclidean lower bound for Pauli coefficients under the
normalized Hilbert-Schmidt norm, applied as a fail-closed feasibility precheck.
It is a conventional baseline optimization. Any reported improvement must
include gate work and state how many orders and grouping proposals were
skipped; it does not establish learning, transfer, or scientific novelty.

## Root implementation audit

`guarded_taylor` retains the same integer derivative and polynomial recurrence.
For `left=T_num^(k+1) epsilon_den` and
`right=T_den^(k+1)(k+1) epsilon_num derivative_den`, the strict tests are
`max|b| > floor(right/left)` and `sum b² > floor(right²/left²)`. These are
exactly equivalent to the rational inequalities. Equality proceeds to the
unchanged proposals. For zero tolerance a nonzero tail is rejected and a
zero tail proceeds to the ordinary zero witness. The denominator of G c_k
is `q0 d^(k+1) k!`, not `(k+1)!`.

Validated numerator entries have at most8192 bits and support at most4096.
The guard rejects unsquared comparison operands wider than16384 bits;
their squares therefore have at most32768 bits, and accumulating at most4096
numerator squares needs at most16396 bits. `max_comparison_bits` records the
unsquared operands. The original exact inputs and bounded order also bound
the power expressions formed before this comparison-budget refusal. The
helper is not a standalone untrusted certificate boundary: the unchanged V7
checker independently validates the accepted output.

The implementation constructs every retained coefficient Fraction before the
gate, but avoids derivative Fraction conversion and all norm proposals for
rejected orders. Early-exit entry inspections/squares are counted. Exact
equality of all accepted polynomials, norm witnesses, bounds and endpoints was
checked across all four arms and all five repetitions in the V11 development
campaign. This supports implementation equivalence in addition to the general
rejection proof; it is not evidence of autonomous acquisition.
