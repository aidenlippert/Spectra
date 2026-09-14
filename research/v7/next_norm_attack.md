# Bounded next norm mechanism: quadratic Pauli-block witnesses

## Scope and conclusion

This is a development-only analysis of the existing `results/v7` examples. I
did not generate held-out cases or change checker code. The usable mechanism is
to certify a small residual sector through its squared Pauli operator, rather
than integrating a coefficient-wise Taylor triangle bound. It can remove a
required Taylor degree on a concrete residual, but the checker work is quadratic
in the retained sector size and therefore is not a claimed headroom win.

The mechanism is distinct from threshold or coefficient pruning: it uses an
algebraic identity and a finite, exact certificate.

## Exact construction

Let a residual sector be

\[
 R(u)=\sum_{i=1}^s q_i(u)P_i,
 \qquad q_i(u)=\sum_{k=0}^d a_{ik}u^k,
\]

with rational Pauli words. The proposer supplies the support and coefficients;
the checker recomputes the residual first. For every ordered pair, it computes
the Pauli product and accumulates the exact coefficient polynomial in

\[
 R(u)^2=\sum_{i,j}q_i(u)q_j(u)P_iP_j.
\]

If all cross terms cancel (for example, a pairwise anticommuting Hermitian
family), the result is the scalar polynomial

\[
 R(u)^2=S(u)I,\qquad S(u)=\sum_i q_i(u)^2.
\]

The checker then proves an outward rational upper envelope `S(u) <= U(u)` on
the interval, using the same exact Bernstein conversion already allowed by the
v7 protocol, and integrates `sqrt(U(u))` by a declared rational envelope. A
simple always-valid version partitions the interval and uses Bernstein
coefficient maxima: if `U_j` bounds `S` on subinterval `j`, charge
`length(j)*sqrt_upper(U_j)`. Square-root intervals are outward and exact.

More generally, a commuting block decomposition gives
`R^2 = diag(S_1 I,...,S_b I)` after the checker verifies the block relations;
the norm bound is the maximum block envelope. Failure of cancellation, block
closure, or an outward square-root enclosure rejects this arm and falls back to
the conventional witness.

## Small exact calculation

Take the two-word sector

\[
 R(u)=(1-u)X+uZ,\quad 0\le u\le1.
\]

Since `XZ+ZX=0` and `X²=Z²=I`,

\[
 R(u)^2=((1-u)^2+u^2)I=(1-2u+2u^2)I.
\]

The coordinate-wise l1 integral is `1 + 1/2 + 1/2 = 2`.
Grouping each power coefficient gives `1 + sqrt(2)/2`. The existing Bernstein
witness already gives 1, since its endpoint coefficients are X and Z. Thus
this example does not establish headroom over the supplied Bernstein baseline. The quadratic witness can use the exact pointwise
norm `sqrt((1-u)^2+u^2)`. A rational partition at `u=1/2` has
`S(u) <= 1` on both halves, hence certifies

\[
 \int_0^1\|R(u)\|_\infty du \le 1.
\]

Thus a tolerance of `1` passes the quadratic witness while the coefficient
triangle witness fails at the same polynomial/order. The bound is conservative
(the exact integral is about `0.812`), and no floating point value is used in
the certificate. In a real residual, this tightening can make order `m` pass
where order `m+1` was previously required; the proposer must report both
attempts and charge the failed attempt under the v7 protocol.

## Resource accounting and comparison

For `s` retained words and polynomial degree `d`, recomputing `R²` costs
`O(s²(d+1)²)` coefficient operations for naive polynomial convolution, plus polynomial envelope work
and square-root enclosure work. The independent checker must repeat the pair
convolution; this is `O(s²)`, not a lower-norm freebie. A block-aware variant
reduces the scalar envelope count only after checking every cross relation.

The conventional alternatives remain coefficient-wise Taylor, BFS projection,
residual-driven sparse expansion, Arnoldi, and existing anticommuting grouping.
The supplied v7 development rows show grouping witnesses often reduce bounds
but add thousands of pair comparisons (for example the `n=4`, XXZ, `gamma=2`,
`T=1/2` order-8 attempt records 2,591 weighted comparisons and still has bound
`0.1231`). Therefore the quadratic witness should be admitted only when its
predicted tightening crosses the tolerance or lowers the required degree by a
full accepted attempt; a lower bound with higher complete construction/checking
cost is not a win.

## Checker schema and refusal rules

The certificate should contain sector labels, exact polynomial coefficients,
the recomputed squared-polynomial coefficients, partition endpoints, and
outward rational square-root intervals. The checker recomputes all of these;
submitted products and bounds are metadata. It must reject non-Hermitian words,
uncancelled cross terms, support leaving the declared invariant block, negative
or malformed rationals, and any interval whose upper endpoint is not proven.
It must charge pair products, coefficient multiply-adds, partition/envelope
operations, rational bit growth, failed attempts, and final checking. This is a
candidate mechanism for bounded replay, not evidence of an acquired method or a
held-out improvement.


## Root development calculation and complete-cost rejection

`experiments/v7_quadratic_probe.py` instead tests actual constant highest-order
residual coefficients of the 24 small-system adaptive Taylor cases. The
quadratic bound tightens 20 and permits one lower Taylor order on two. The
four-qubit case requires 32,385 unordered pair tests for the square; norm
improvement alone is not economical evidence.

A separate exact checker, corruption tests, and seven randomized paired repeats
then test selective squaring across the 36 unchanged development problems.
The only order-saving eligible case is the three-qubit XXZ case at gamma=2,
T=1/2: degree 12 drops to 11, but complete median runtime increases by about
11.3%. An additional eligible case pays for a failed square attempt. This
candidate is rejected as an efficiency mechanism in its current implementation.
The original residual checker and held-out distribution remain unchanged.
