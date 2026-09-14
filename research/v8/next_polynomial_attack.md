# Next polynomial attack: rational Chebyshev collocation (preflight only)

The segmented-grid candidate is closed: reducing the interval multiplies exact
support generation and introduces charged endpoint jumps.  The next distinct
question is whether a single non-Taylor polynomial can lower total cost under
the same `v7-residual-1` checker.

## Precise candidate

Use degree `m` Chebyshev interpolation of the observable trajectory on one
interval.  Let `x=2t/T-1`, choose the rational Chebyshev nodes

`x_j = cos((2j+1)pi/(2(m+1)))`,

and replace these irrational nodes by an exact rational node set, for example
`x_j = -1 + 2j/m` (the resulting basis is still a Chebyshev polynomial basis,
but the interpolation grid is rational).  At each node construct the value
`Y_j` by a bounded numerical Krylov/exponential proposal, then rationalize all
coefficients by an outward interval enclosure.  Solve the Vandermonde system
exactly (fraction-free Bareiss elimination) to obtain Chebyshev coefficients
`a_0,...,a_m`; convert the result exactly to the power coefficients expected by
`Piece`.  The exported candidate is only this rational power polynomial.  The
checker does not trust node values, Krylov values, or rounding.

This precise choice is intentionally conservative: a genuine Chebyshev node
implementation would require algebraic-number certificates, so rational nodes
avoid an unimplemented new checker.  It is also ordinary polynomial
interpolation, a known numerical method rather than a proposed learner.

## Certification lemma and cost obligations

For `P(t)=sum_k c_k t^k`, define `R=P'-GP`.  The existing checker certifies
`||P(0)-O|| + integral_0^T ||R(t)|| dt` by exact residual records.  In power
coordinates it uses the monomial integral; in Bernstein mode it forms the
exact degree-normalized combinations already implemented in
`residual_records`.  Therefore a Chebyshev proposal needs no new functional
checker: after exact basis conversion, the unchanged checker proves the same
all-state final-time observable bound.  The interpolation theorem is not a
certificate; it can guide node selection only.

The full cost must include, for every node, construction of `Y_j` and its
independent numerical enclosure, the dense `(m+1)x(m+1)` fraction-free solve,
Chebyshev-to-power conversion (coefficient growth can be exponential in `m`),
all generator applications during checker replay, residual witness grouping,
and failed orders.  A hidden full-matrix exponential or an uncharged solve is
invalid.  With `m<=24`, Bareiss arithmetic is bounded but may dominate sparse
Taylor: Taylor needs `m` sparse recurrence applications, whereas collocation
needs `m+1` trajectory solves plus `O(m^3)` exact elimination and then still
requires the checker to recompute `m+1` residual coefficients.

## Small algebraic preflight

The two-dimensional rotation generator `G(X)=-Y, G(Y)=X`, initial observable
`X`, is an exact sanity case.  A degree-one interpolation using the first two
terms of the node-value proposal gives `P=X-tY`; its residual is directly
recomputed and is not exact dynamics.  Degree-two
collocation at `t=0,T/2,T` gives a rational polynomial and can be compared
against Taylor degree two using only `Fraction`, `Generator.apply`,
`derive_certificate`, and `check_certificate`.  This preflight is bounded and
meaningful because it exercises basis conversion, exact solve, endpoint jump,
and the full-horizon residual bound.  It must report both accepted bound and
all solve/check work; a smaller residual alone is insufficient.

For this rotation, Taylor already has the exact recurrence support `{X,Y}` and
zero checker support growth.  Collocation cannot reduce generator work and
adds the dense solve, so it cannot be a credible win on this sanity case.  A
nontrivial 2/3-qubit preflight should therefore be attempted only if the
small case shows no arithmetic pathology; it cannot be used to claim a win on
the protocol distribution.

## Conditional headroom decision

The only plausible win is a case where endpoint/node trajectory values have
much smaller sparse support than all Taylor powers and the exact solve remains
cheap.  That condition is not established by the current artifacts.  For the
finite-spin generators, nonnormal commutator dynamics also means a scalar
Chebyshev spectral interval does not provide a sound error bound; the residual
checker remains essential and removes the usual black-box Chebyshev remainder
shortcut.  Arnoldi already supplies a stronger conventional projected
polynomial baseline, including explicit orthogonalization and solve counters.

Consequently this is a conditional research direction, not a promised
headroom result.  Any implementation must first pass the bounded 2/3-qubit
preflight and then a full common-case comparison against fraction-free Taylor
and actual v7 Arnoldi, charging every solve, rationalization, residual replay,
and refusal.  No held-out evaluation or method learner is authorized by this
preflight.
