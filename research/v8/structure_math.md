# V8 structure attack: symplectic parity/block recurrence (negative)

## Candidate representation

Encode a Pauli word by the binary symplectic vector `v=(x|z)`, with `Y=(1,1)`.
For a Hamiltonian word `q` and current word `p`, a nonzero commutator is exactly
the edge `p -> p xor q` when `<p,q>=1`, with the sign determined by the
symplectic phase.  Thus the sparse generator is a weighted Cayley graph on
`F_2^(2n)`.  A possible exact acceleration is a linear charge `l(v)` such
that `l(q)=1` for every Hamiltonian term q.  The commutator then toggles charge
at every edge, and a Taylor recurrence can be represented as alternating
even/odd blocks (or, more generally, as a block recurrence for a small set of
charge sectors).  This is an algebraic representation, rather than cache or
gcd reuse.

## Development probe

`results/v8/structure_probe.json` enumerates all `2^(2n)` linear charges for
the prescribed Hamiltonian terms, and separately explores the actual
commutator graph from the central-Z seed.  No global linear charge toggles on
all Hamiltonian terms.  However, the *reachable* graphs at n=3 and n=4 are
fully explored and bipartite (30/63 and 126/255 vertices for XXZ/mixed).  The
bipartition is explained exactly by a universal real/imaginary structure:
every workload Hamiltonian word has even Y-parity, hence H is real.  A
Hermitian Pauli is real or purely imaginary according as
`c(P)=#Y(P) mod 2`, and `i[H,P]` swaps these sectors.  Thus
`c(x|z)=sum_i x_i z_i mod 2` is the nonlinear symplectic charge.  n=6 hits
the 512 support cap before closure and is explicitly marked incomplete; its
reported bipartite status is therefore unknown, not a conclusion.

## Why the general representation does not create headroom

Without a uniform charge, one can still compute the complete symplectic edge
test and group words by the values of several linear charges.  This requires
the same Hamiltonian-term/Pauli-pair commutation tests already performed by
`Generator.apply`; each accepted edge still needs a coefficient multiply-add.
This real/imaginary bipartition is a known time-reversal structure and is
already implicit in the Pauli labels and commutator routine.  Applying it
requires only a Y-parity bit per word, but does not remove any nonzero
commutator edge or coefficient multiply-add.  Damping is diagonal and
preserves the sectors, so the full generator has off-diagonal parity blocks
plus same-sector diagonal damping.  A block norm identity could be explored,
but would require a new checker argument and is not a free construction saving.
Building a sector index costs `O(|H| |B|)` binary dot products and stores an
extra sector key per live word.  Applying a block recurrence costs the same
number of nonzero weighted edges as the ordinary sparse recurrence, with
additional sector dispatch and, for two-step blocks, intermediate support
materialization.  It cannot beat the strongest conventional baselines, which
already retain all generated support, cache columns, and use exact rational
recurrences.  A special case where all terms share a charge is not present in
the frozen development families and would be an instance-specific symmetry,
not a reusable learned method.

An exact two-step recurrence `G^2 b_k` is also not cheaper: it performs all
edges leaving the first frontier and then all edges leaving the second; unless
many paths cancel, it strictly adds intermediate coefficient work.  Detecting
and proving cancellation requires exact accumulation of those paths, so it
does not remove the baseline work.  If a block is formed explicitly, the
checker must still recompute `G c_k` for every exported coefficient and verify
all residuals, so the representation cannot reduce independent checking cost.

## Proof obligations for any future symmetry claim

1. Prove the charge equations over every supplied Hamiltonian term and include
   zero/degenerate sectors; no numerical inference is sufficient.
2. Prove that each block recurrence is algebraically identical to `i[H,·]`
   with exact signs and damping rates, including collisions where multiple
   terms map to one Pauli word.
3. Prove support and rational-bit budgets, and charge index construction,
   collision checks, sector conversion, failed attempts, and validation.
4. Export the ordinary coefficient maps and pass the unchanged independent
   residual checker; a sector transcript cannot be trusted by the checker.
5. Compare complete construction + witness + checker cost against full Taylor,
   residual-adaptive, BFS projection, and Arnoldi baselines.  A lower edge
   count on a retained representation alone is insufficient.

## Decision

No cost-positive reusable mechanism was demonstrated.  The parity structure is
an exact explanation of the observed graph but supplies no demonstrated
operation-count reduction against existing sparse baselines.  No learner or
held-out n=5/7 generation occurred, and this report makes no conditional
theorem or learned-capability claim.
