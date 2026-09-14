# V10 dense matrix baseline: cost and certification audit

## Finding

A dense baseline is a legitimate missing conventional comparator for the small
systems.  The smallest credible version is a matrix-free dense `d x d`
representation of the *state-independent observable generator*, followed by
a certified outward-rounded action on the initial operator.  It must be
treated as a bounded comparison for small `d`, not as a general solver.

The current timing motivates a fresh paired measurement: reference Duhamel is
about `0.241 s`, while fraction-free Taylor is about `0.0159 s` on n3 XXZ,
`gamma=2`, `T=.5`.  A dense arm could still be useful as a conventional
cross-check or for a very small case; its cost and certification must be
measured contemporaneously with the matched Taylor call.

## Smallest algorithm

Let `d=2^n`.  Do not construct a `q x q` Liouville matrix.  Store the exact
observable as a dense `d x d` matrix and apply the state-independent generator
matrix-free:

`G(O)=i(HO-OH) + (gamma/4) sum_{i,a}(sigma_{i,a} O sigma_{i,a}-O)`.

Here the Pauli conjugations are signed permutations, so each damping term costs
`O(d^2)` rather than a generic matrix multiplication.  A dense Hamiltonian
left/right multiplication costs `O(d^3)` per action (or `O(L d^2)` when the
supplied `L` Pauli terms are applied as sparse signed-permutation maps).  For a
dense comparison, represent `H`, `O`, and the recurrence arrays as integer or
dyadic matrices and compute either:

* a truncated Taylor action `sum_{k=0}^m T^k G^k(O)/k!`, with an exact rational or
  dyadic residual bound; or
* a dense `expm`-style action only if its backward/rounding certificate is
  supplied.  A black-box floating eigensolve is a timing cross-check, not a
  certificate.

The first option is the smallest sound baseline because the existing proof
already gives, under the declared contractive semigroup assumption,

`||U(T)O-P_m(T)||_infty <= ||O||_infty (K T)^(m+1)/(m+1)!`,

with `K >= ||G||_{infty->infty}` (or a separately proved tighter bound).  The
dense matrix is then only an implementation of the same certified action;
it does not get to replace the physical contractivity assumption with a
numerical spectral-radius estimate.

For the exact initial operator, form `O` from rational Pauli coefficients and
Kronecker products of `I,X,Y,Z`; count this conversion.  If the final result
must be reported in Pauli coordinates, use the exact trace projection and
charge it too.  The matrix-free route has `d^2` storage and avoids any
`q^2` Liouville allocation.

## Complete cost accounting

The charged calculation must include all of the following, including failed
orders and refusals:

1. Hamiltonian and initial-matrix construction, including exact Kronecker
   products. Applying one generator action costs `O(d^3+n d^2)` for dense `H`,
   or `O(L d^2+n d^2)` for sparse Pauli terms.
2. Each Taylor step applies `G` to a `d x d` matrix, so degree `m` costs
   `O(m(d^3+n d^2))` (or the sparse bound above). The denominator/factorial
   and residual arithmetic are included.
3. Range checks and arithmetic-error bounds for every multiply/add. The
   peak absolute entry bound `B_j` must be recorded per step; for a dense
   matrix product a safe entrywise recurrence uses `d B_A B_H` (and the
   analogous right-product bound), plus `n` signed-permutation damping terms.
   Overflow or a failed bound is a refusal.
4. Certified endpoint and readout. If the output remains in Pauli
   coordinates, the observable error is bounded directly by its coefficient
   l1 residual. If a computational-basis matrix is used, Pauli projection,
   Hermiticity checks, and the matrix-to-Pauli l1 readout are charged.
5. Independent checker replay, serialization, and all rejected candidate
   orders, with timing separated only for reporting—not removed from complete
   cost.

The ambient Pauli counts remain `q=63,255,4095` for n=3,4,6, but they are
readout/support counts only.  The matrix-free physical matrices have
`d=8,16,64` and storage `d^2=64,256,4096`; n7 has `d=128` and `16384`
entries.  This does not make the method universally cheap: the cubic dense
Hamiltonian action and arithmetic precision still grow, but no `q^2`
Liouville matrix is formed.

## Can NumPy integer/dyadic arithmetic be rigorous?

Not by itself. NumPy fixed-width integer multiplication is rigorous only when
every intermediate is proven to stay inside the selected dtype range. Matrix
multiplication can overflow silently, and a post hoc cast cannot recover the
lost high bits. `dtype=object` with Python integers avoids fixed-width overflow
but gives Python big-integer arithmetic and does not establish a faster
baseline than `Fraction` without measurement.

A rigorous faster variant is possible only with a fixed dyadic format and an
explicit proof before each kernel that all exact numerators fit the chosen
width, together with outward rounding after every operation.  For an `s`-bit
signed integer matrix, if `A` and `B` have entry magnitudes bounded by `M_A`
and `M_B`, a matrix product of inner dimension `d` is bounded by
`d M_A M_B`; an accumulated sum of `r` such products is bounded by
`r d M_A M_B`.  For the sparse Pauli action replace `d` by the number of
nonzero signed-permutation contributions.  The Taylor recurrence must carry
these bounds step by step, including division by `k` and the chosen dyadic
rounding radius.  Every product and sum gets an outward error interval; the
next-step radius is bounded by the same absolute-value action norm applied to
the prior radius plus the fresh kernel roundoff.  Reject before a kernel if
the proven numerator or radius can exceed the dtype range.  `dtype=int64`
without this forward bound is not rigorous.  Arbitrary-precision integers
with `[lo,hi]` dyadics are rigorous, but any speed advantage over `Fraction`
remains empirical.

Exponentials are not needed for the Taylor baseline. If an exponential or
trigonometric endpoint is used, it requires a separately verified outward
dyadic enclosure (range reduction plus a Taylor/continued-fraction bound).
Calling NumPy floating `exp` and widening by a guessed epsilon is not a
certificate.

## Error budget

Use the fixed total tolerance `1/1000`. Allocate, for example,

`epsilon_tail + epsilon_arith + epsilon_readout <= 1/1000`,

with each allocation represented by exact dyadics. `epsilon_tail` is the
factorial bound above. `epsilon_arith` is the accumulated interval radius of
the dense recurrence, and `epsilon_readout` covers Pauli projection or
observable evaluation. The initial operator is exact, so it contributes zero
initialization error; if it is converted through dyadic matrices, the exact
conversion must be checked and any nonzero rounding charged. Hermiticity and
the damping diagonal must be checked independently. A numerical residual
alone does not certify the full semigroup action.

## One next preflight

Run exactly one paired case: n3 XXZ, `gamma=2`, `T=.5`, central-Z observable,
`epsilon=1/1000`. Apply `G` to checked dense physical matrices without building
the Liouville matrix;
try the first factorial-valid orders up to 24; count construction, every
failed order, interval arithmetic, endpoint readout, and independent replay.
Require a certified total error at most `1/1000`, and compare complete costs
in the same contemporaneous paired run against fraction-free Taylor.  Include
construction, every failed order, interval arithmetic, endpoint readout, and
independent replay in both arms.  A pass on n3 is only a reason to repeat on
the already admitted n4 and n6 cases; do not infer n5/n7 performance or a
universal solver from it.
