Exploratory conditional route; no new factorization was constructed. See DERIVATION.md for the implemented spectral witness.

# Compact closure routes for the positive-cone certificate

The current amplitude certificate is exact but stores a dense `70\times70`
matrix `C` and a dense residual `D=L(C)-\ell C`.  A compact `C` representation
alone is not enough: the accepting condition is the global PSD inequality
`D\succeq0`.

## A usable tensor-cone theorem

Let the amplitude space be a one-dimensional chain of local factors. Suppose
there is an exact factorization

\[
 C=VV^*,
\]

where `V` is given as an MPS/MPO with bond dimension `r`, and suppose the
residual admits an exact local Gram factorization

\[
 D=L(VV^*)-\ell VV^*=\sum_{a=1}^k R_a^*R_a,
\]

with each `R_a` represented by an MPO of bond dimension at most `s`. Then
`C\succeq0` and `D\succeq0` follow without forming either dense matrix. Exact
replay only needs to verify the tensor contraction identities and the local
MPO coefficient identities. If `V` has full row rank, then `C\succ0` follows
from a certified lower singular-value bound for `V`; a numerical rank estimate
is not sufficient.

This is a valid compact positivity route, but the hard condition is the second
factorization. Applying a local map `L` to a low-bond `V` can increase bond
dimension, and regrouping the result into `R_a^*R_a` can be as hard as the
original PSD problem.

## Boundary-contraction variant

For a chain ordered by cuts, a stronger acceptance certificate is a recursive
Schur complement: expose one tensor layer at a time, certify each local Schur
complement as PSD, and carry a positive boundary matrix of dimension `b`. If
every boundary matrix has an exact rational PSD factorization and
`b` remains bounded, induction proves global PSD. The proof cost is then
`O(n poly(b,s))` tensor contractions rather than dense dimension squared.

The induction must include every boundary cross term. Dropping a small tensor
coefficient without a norm allowance does not preserve positivity. A practical
checker should replay (i) local coefficient equality, (ii) boundary update
identity, (iii) PSD of each local Schur complement, and (iv) a final boundary
PSD condition.

## Why ordinary tensor compression is insufficient

Low bond dimension of `V`, short correlation length, or an area law does not
imply that `D` has a low-bond Gram factorization. Positivity is a cone
constraint, not an approximation norm. A low-bond matrix can have a local
operator image whose positive and negative terms cancel globally, requiring a
large bond dimension to expose a sum of squares. Higher-dimensional tensor
network contraction is also hard in the worst case, so the claim must be tied
to the actual one-dimensional amplitude ordering and an explicit boundary
bound.

Likewise, a purification or MPDO with small bond dimension certifies `C\succeq0`
but says nothing about `L(C)-\ell C\succeq0` until the residual cone is exposed.

## Acceptance theorem for a proposed compact certificate

A compact candidate should be accepted only if it supplies exact rational
local tensors and a replayable proof of:

1. `C=VV^*` and `C\succ0` (or an explicit positive lower singular-value bound);
2. `D=L(C)-\ell C` coefficient equality after all fermionic signs and
   chemical-shift identities;
3. a local Gram or recursive Schur certificate for `D\succeq0`;
4. a bound on tensor bond dimensions, boundary dimensions, and cold
   construction/replay costs.

Post-hoc SVD compression of dense `C` or `D` is a diagnostic only. It does not
establish cheap discovery or exact positivity. The strongest near-term test is
to seek such a factorization for the verified H8 candidate, then repeat it on
two coupled blocks; failure should be reported as failure of the selected
tensor cone, not as an impossibility theorem for positive-cone methods.
