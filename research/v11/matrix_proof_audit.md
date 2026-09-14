# V11 matrix proof audit

This is a bounded proof review of the proposed exact matrix-action route. It
does not establish novelty or any performance/scientific-transfer claim. The
comparison target is the existing V8 Taylor polynomial and the V7 checker.

## 1. Exact operator and denominator clearing

For `d=2^n`, retain `O in C^(d x d)` and apply

```text
G(O) = i[H,O] + gamma/4 * sum_(j=1..n,a in {X,Y,Z})
       (sigma_(j,a) O sigma_(j,a) - O).
```

On a Pauli string `P`, each one-qubit conjugation is either `P` or `-P`.
Therefore the damping term is
`(gamma/4) * sum_(j,a) (sign_(j,a,P)-1) P = -gamma*weight(P)*P`:
for each nonidentity letter exactly two of the three axes anticommute, giving
`-2*(gamma/4)` per occupied qubit.

Choose a positive integer `q` clearing every coefficient of `H` and
`gamma/4`. Write `h_r=q*H_r` for the cleared Hamiltonian coefficients and
`g=q*gamma/4`. Then `D=qG` has Gaussian-integer matrix entries whenever the
input matrix has Gaussian-integer entries. This `q` is a coefficient
denominator; it is not a Liouville-space dimension.

## 2. Signed-permutation action and component bounds

For every Pauli word `R`, left multiplication by `R` maps a computational
basis matrix entry to one entry with phase in `{1,-1,i,-i}`; right
multiplication does the same. Hence `R O` and `O R`, and consequently
`R O R`, are signed/permuted copies of the real and imaginary integer arrays
(with a possible real/imaginary swap). A Hamiltonian commutator term has two
such copies. If every real and imaginary component of `O` is bounded by `M`,
the componentwise bound for one application of `D` is

```text
M_next <= (2 * sum_r |h_r| + 6*n*|g|) * M.                 (1)
```

The damping count is `3n` conjugations plus `3n` subtractions, each weighted
by `|g|`; (1) intentionally uses the safe `6n|g|` bound. It bounds temporary
accumulators as well as the final cancelled output only if the implementation
checks the running sum after every add. The same bound can be iterated for
each derivative, with a separate bound for the current order's coefficient.

Before entering an int64 kernel, use Python/unbounded integers to check all
of the following:

* every input component and every scaled coefficient;
* every product `component * |h_r|` and `component * |g|`;
* every per-entry Hamiltonian/damping accumulator, using (1) (or the exact
  term-by-term sum when tighter);
* every derivative order and Taylor numerator, including factorial and `q`
  scaling;
* every Walsh butterfly row sum, bounded by `d*M` for that row.

The check must use the largest temporary magnitude, not the final value after
cancellation. If any checked bound is outside signed int64, dispatch to an
explicitly costed arbitrary-integer path or refuse. Detecting overflow after
fixed-width arithmetic is invalid. Track peak temporary magnitude, final
matrix magnitude, rational bit length, matrix-entry count, and Pauli-support
cap under the same budget policy as V7/V8.

## 3. Fast Pauli conversion

Use bit masks `x,z` and
`P(x,z)=i^popcount(x&z) X^x Z^z`, with the local convention `Y=iXZ`.
For computational basis state `|s>`,

```text
P(x,z)|s> = i^popcount(x&z) (-1)^popcount(z&s) |s xor x>.
```

Expanding `trace(P A)` therefore gives the exact identity

```text
c(x,z) = trace(P(x,z) A)/d
       = i^popcount(x&z)/d
         * sum_s (-1)^popcount(z&s) A[s, s xor x].       (2)
```

For fixed `x`, the inner sums in (2), over all `z`, are one Walsh-Hadamard
transform of the diagonal-offset row `A[s,s xor x]`. The transform uses only
integer additions/subtractions on the real and imaginary arrays; apply the
phase in `{1,-1,i,-i}` and divide by `d` exactly. Its row sum is bounded by
`d*M` before cancellation. This is a standard transform and a rederivation
for this audit, not an autonomous-discovery or novelty claim.

## 4. Equality with the V8 Taylor polynomial

Let `q0` clear the initial matrix/Pauli coefficients and let `b_0=q0 O_0`
be the integer initial numerator. Let `b_{k+1}=D(b_k)` where `D=qG`. The
matrix coefficient exported at order `k` is

```text
c_k = b_k / (q0 * q^k * k!).
```

Applying `G` to it gives the derivative coefficient

```text
G(c_k) = b_(k+1) / (q0 * q^(k+1) * k!).                 (3)
```

Equation (3) is the required derivative coefficient; the `(k+1)` in the
polynomial derivative is supplied separately by `(k+1)c_(k+1)`. Thus the
matrix recurrence and the V8 integer recurrence produce the identical rational
Taylor polynomial after conversion by (2). Their integer clearing factors can
differ: the matrix code clears gamma/4, whereas V8 clears gamma. For gamma1/5
the matrix factor can be20 while V8 uses10. Each numerator carries its own
matching denominator; equality is asserted only after exact rational conversion.
This equivalence requires the same
Hamiltonian convention, damping `gamma/4`, Pauli/Y phase convention, term
ordering only for cost accounting, and exact arithmetic. It must be checked
by exact coefficient equality on bounded n1..4 cases before timing claims.

## 5. Checker and cap semantics

After each exact derivative, convert to Pauli coefficients, construct the same
power-basis tail witness, and charge conversion work. Reuse V7's immutable
checking rules unchanged: exact rational parsing and bit caps, term/support
caps, required record coverage, duplicate-free term coverage, pairwise
anticommutation checks, exact square sums, nonnegative upper bounds, exact
claimed-bound equality, and requested-time equality. A conversion or action
budget failure is a rejection/refusal, never a partial certificate.

The implementation must preserve the same `max_terms`, rational bit-length,
matrix-entry, and order caps for both the cached checker path and the live
original path. “Checker512” cached and live-original runs must be compared
under identical cap semantics and accounting: cache hits may reduce repeated
Hamiltonian construction cost, but cannot waive action, conversion, norm, or
pair-check charges. Exact norm checks must recompute the supplied witness
against converted coefficients; a faster cached result is not evidence that
the original path accepted a different object.

## 6. Narrow acceptance conditions

The first paired case remains V10's existing `n=3` XXZ, `gamma=2`, `T=1/2`,
`epsilon=1/1000`, central-Z observable. Acceptance of this audit requires:

1. exact coefficient equality with V8 on bounded small cases;
2. pre-kernel integer bound checks covering products, partial sums, factorials,
   Walsh rows, and temporary as well as final magnitudes;
3. explicit arbitrary-integer fallback or closed refusal when int64 is unsafe;
4. identical V7 checker semantics and caps on cached/live-original paths; and
5. all endpoint conversion and fresh exact checking included in the cost.

If the complete method loses the paired cost comparison, retain the exact
conventional result and identify whether action, conversion, or immutable
checking dominates. A win licenses only development expansion and paired
repetitions; it does not open heldout evaluation or establish learning.
