# V9 norm subdivision: a bounded certificate idea and its cost obstruction

This note evaluates one narrowly scoped extension of the V7/V8 residual
certificate.  It is a mathematical preflight only; it is not an acquisition
result or a claim about the intended physical system.

## Existing baseline

`residual_records` recomputes the exact sparse residual coefficients.  In
power mode it charges

\[
  \sum_{k=0}^m {T^{k+1}\over k+1}\,\|r_k\|,
\]

and in Bernstein mode it applies the exact power-to-Bernstein transform and
charges `T/(m+1)` times the norm of every transformed coefficient.  Thus the
Bernstein option already has the right integral weights, but it still applies
the triangle inequality coefficient by coefficient.  It cannot see a sign
change or cancellation between Bernstein coefficients.  The V8 development
data show the practical scale: a 3-qubit, degree-5 candidate has 53
polynomial entries and about 0.002 s checking time, while a 4-qubit,
degree-13 candidate has 1,175 entries.  A separate profiling artifact reports
larger replay timings, but that includes profiling overhead and is not the
ordinary V8 complete-check time.  Width 6 already reaches the 512/1024
cached-column caps in several cases.

## Exact subdivision bound

Let an operator-valued residual on one piece be represented in the degree-
`m` Bernstein basis on `[0,T]` by `C_0,...,C_m`.  For a dyadic interval `I`
of length `h`, run de Casteljau at the rational endpoints of `I` and retain
the local Bernstein control operators `C_{I,0},...,C_{I,m}`.  Since the
Bernstein basis is nonnegative and sums to one,

\[
  \|R(t)\|_\infty \leq \max_i\|C_{I,i}\|_\infty\quad(t\in I).
\]

Consequently the finite, independently checkable bound

\[
  \int_0^T\|R(t)\|_\infty dt
  \leq \sum_{I\in\mathcal P} |I|\,
       \max_i U(C_{I,i})
\]

is sound whenever each `U` is an independently checked Pauli norm upper
bound.  A checker need only verify the de Casteljau affine recurrences over
the rationals, verify each local norm witness, and compare the advertised
sum.  No quadrature and no floating point estimate are involved.  A slightly
stronger, still simple variant uses the coefficientwise integrated envelope
`sum_i |I| U(C_{I,i})/(m+1)`, which can be smaller than the interval maximum.

The useful theorem is therefore: for every finite subdivision `P`, the above
quantity is an upper bound.  At a split `u`, de Casteljau uses
`d_i^(r+1)=(1-u)d_i^r+u d_(i+1)^r`, so child controls are convex
combinations of parent controls.  Norm convexity plus preservation of the
Bernstein coefficient mass proves that exact-norm integrated envelopes cannot
increase under refinement.  With independently rounded upper witnesses,
retain the parent candidate and take the minimum, since arbitrary rounding
need not preserve monotonicity.  The theorem is
purely a Bernstein convex-hull statement and does not use diagonalization,
normality, or a matrix norm estimate hidden in a numerical routine.

## Finite cost and why it is unlikely to buy degree headroom by itself

With degree `m`, `s` subdivision levels, and `N` distinct Pauli labels per
control operator, a straightforward exact checker performs

\[
  O(N m 2^s)\quad\text{stored coefficient visits},\qquad
  O(N m^2 2^s)\quad\text{rational multiply/add work}
\]

if every child row is formed independently.  Shared de Casteljau rows reduce
the constant but not the `2^s` number of interval witnesses.  It also needs
`(m+1)2^s` norm witnesses per residual polynomial (or a transcript that
reconstructs them), versus `m+1` in current Bernstein mode.  At degree 24,
four levels already mean 400 local controls per residual, before Pauli-group
pair checks and exact rational bit growth.  Reusing the original generator
cache helps only the initial residual construction; it does not remove the
new norm and scalar-bound work.  The independent checker must still replay
all residual coefficients, so subdivision cannot erase the dominant replay
cost documented for V8.

The method is attractive only when a very small number of subdivisions gives
a large bound reduction.  A deterministic admission rule should cap `s` and
charge every generated child, rational operation, norm witness, and checker
comparison.  It should also reject a subdivision certificate whose claimed
bound is not exactly the sum of the checked interval terms.

## Explicit cancellation example

For the scalar residual `q(t)=1-2t` on `[0,1]`, the degree-one Bernstein
controls are `(1,-1)`.  Current coefficientwise Bernstein integration gives
`(1/2)(|1|+|-1|)=1`, although the exact integral of `|q|` is `1/2`.
Subdividing at `1/2` gives child controls `(1,0)` and `(0,-1)`.  The
coefficientwise integrated envelope is
`(1/2)[(1/2)(1+0)] + (1/2)[(1/2)(0+1)] = 1/2`, exactly sharp.  This is realizable
as a single residual coordinate of a scalar rotation calculation: take a
Pauli `X` coefficient `a(t)=t-t^2` in a candidate, so its derivative residual
coordinate is `a'(t)=1-2t`; for `H=Z/2` the remaining commutator coordinate is
checked separately.  The example proves finite cancellation headroom in one
coordinate, but it does not show a full operator certificate below `10^-3`.

## Decision for V9

Subdivision is mathematically valid and could be a useful bounded diagnostic
for low-degree, low-support residuals.  Existing evidence does not justify
assuming it will rescue degree `<=24`, widths 3/4/6, or tolerance `10^-3`:
the scalar gain requires a sign-changing residual, while Pauli norm witnesses
sum distinct operator directions and can leave no cancellation to expose.
The appropriate next test is a development-only replay on the existing
degree/order grid, with `s <= 2` first and complete construction, witness,
subdivision, and independent-check costs recorded.  Unless that replay shows
a certified case whose total charged cost beats the current Taylor/Bernstein
baseline, subdivision should remain a bounded mathematical extension rather
than a promoted construction rule.
