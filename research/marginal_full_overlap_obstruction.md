# Exact full-overlap obstructions beyond the three-spectator family

The selected ENERGYv16/FAMILYv12 local mixtures fail complete five-site overlap
consistency in both W=0 and W=1 cases. Exact positive-projector probabilities
prove this, without numerical eigenvalue acceptance. Each witness also gives
a compact, symmetry-compatible stationary telescope. The two telescopes are
linearly independent modulo the preceding affine operator family.

These findings refute a translation-invariant quantum extension of these
particular symmetry-averaged six-site mixtures. They do not invalidate the
accepted family ceilings or energy lower bounds. No new energy certificate or
family ceiling was produced in this diagnostic turn.

## Complete reduced density matrices

Inputs are the 137-source positive, normalized mixtures in
`results/marginal_graded_hubbard8/three_spectator/<case>/final/`.
Each integer source is averaged over its eight particle-hole, spin-flip and
reflection images. The implementation checks closure of each projective orbit
and reconstructs both contiguous five-site partial traces with a common exact
integer denominator. Unit trace and spin-sector support are checked exactly.
Matrix hashes, exact rational quantities and source hashes are retained in
`results/marginal_graded_hubbard8/full_overlap/<case>/full_overlap_replay.json`.

Let rho_L and rho_R be the left and right reduced states. For a vector v, use
the normalized positive projector P=vv*/||v||². A stationary extension requires
Tr(P rho_L)=Tr(P rho_R). The following displayed decimals summarize exact
rational results in the receipts.

| Quantity | W=0 | W=1 |
|---|---:|---:|
| Witness v | \|358> - \|601> | \|346> + \|613> |
| Squared vector norm | 2 | 2 |
| Five-site spin numbers | (3,2) | (3,2) |
| Left probability | 0.03455734832389095 | 0.00782211122115541 |
| Right probability | 0.034272965597014836 | 0.00807266225122646 |
| Left minus right | +0.0002843827268761191 | -0.0002505510300710504 |
| Nonzero upper-triangle differences | 30,952 | 31,000 |
| Nonzero diagonal differences | 648 | 664 |
| Squared Frobenius norm of difference | 0.000008287674679483734 | 0.000020589858053225703 |
| Common denominator bits | 4,709 | 4,926 |

Sites are encoded with two fermionic bits each, with site zero in the low bits.
The four five-site basis states, in site order, are:

- 358: down, up, down, up, up.
- 601: up, down, up, up, down.
- 346: down, down, up, up, up.
- 613: up, up, down, up, down.

Every site is singly occupied. Each witness pair differs in eight fermionic
bits: two simultaneous spin exchanges. The diagonal differences concern full
spin-resolved probabilities. They do not contradict the previously verified
stationarity of the coarser charge-only laws or their classical Markov extensions.

All four probabilities lie in [0,1]. The violation is unequal overlapping
probabilities, not a negative probability or a failure of local positivity.
The absolute mismatch also lower-bounds the trace distance between the two
five-site marginals because P is an effect between zero and the identity.

## Compact stationary telescopes

Let H, S and R denote the signed five-site particle-hole, spin-flip and
reflection transformations. With P as above, define

    A = (P + H P H* + S P S* + S H P H* S*) / 4
    Y = (A - R A R*) / 2
    T = Y_left - Y_right.

Five-site particle-hole squares to minus the identity on vectors; its density
action is projectively well-defined. This sign is explicitly tested rather
than assuming a vector involution. The numerator of Y has only 16 nonzero
entries in either case, with common denominator 16. The embedded T has 116
nonzero entries for W=0 and 124 for W=1 in the full 4096-state Fock space.

Exact checks establish Hermiticity, spin-number conservation and invariance of
T under all three six-site symmetries. Fermionic cyclic translation of the left
embedding equals the right embedding, including fermionic permutation signs.
The six translated copies cancel exactly on the full six-site Fock space.
The local embedding identity gives the same telescoping cancellation on larger
periodic chains; on an open chain only the two endpoint embeddings remain.

Since A and RAR* are positive contractions, ||Y|| <= 1/2 and ||T|| <= 1.
The open-chain endpoint difference likewise has norm at most one. These are
analytic norm bounds from normalized projectors; numerical eigenvalues are
used only in an independent test. Direct exact expectation of T on the original
137-source mixture agrees with the full symmetry-averaged RDM mismatch above.

The receipts are `<case>/telescope_replay.json`. The intermediate
`projector_witness.json` files remain marked nonaccepting proposals; the later
telescope receipts accept only the stated algebraic and consistency claims.

## Independence from the previous operator family

The new W=0 telescope has matrix element -1/8 at (1382,1625); the W=1 telescope
has +1/8 at (1370,1637). Both pairs lie in the six-site spin sector (4,2) and
change eight fermionic bits. The two new operators restricted to these entries
have a two-by-two determinant of -1/64.

An independent exact support replay checks all 71 preceding nondiagonal
telescopes on all 4096 Fock states: one hopping, four spin, fourteen
one-spectator hopping, four pair-transfer, thirty two-spectator hopping and
eighteen three-spectator hopping terms. They change at most four bits, and all
vanish at the two selected entries. The actual half-filled projector sources
are supported in (3,3); the charged projectors have total particle number five
or seven. They also vanish at these (4,2) entries. Diagonal terms and physical
one-body hopping cannot contribute there.

Thus the two telescopes add two independent directions modulo the full prior
affine family. This is an operator-span statement, not strict separation of
energy bounds or a claim that adding these two conditions closes full overlap
consistency. The supporting receipt is `full_overlap/support_replay.json`.

## Validation and remaining questions

All 21 focused tests pass, covering independent dense partial traces and
Kronecker embeddings, fermionic translation via creation operators, projective
particle-hole signs, direct versus reduced-state moments, valid probabilities,
stationary examples, rescaling and invalid source refusals. Both telescope
replays and the prior-support replay run with `python -S` using only the
standard library. The full-RDM acceptance code likewise uses only standard
library integer and rational arithmetic.

No production certificate implementation changed this turn. The earlier full
suite result remains 982 tests plus 102 subtests, with three subsequent
fraction-free discovery tests run separately. This turn's 21 focused tests
are not a new full-suite run. The collector audits the earlier 479-file
provenance manifest and all five new accepted receipts against current hashes.
All launched jobs are terminal; no GPU, paid resources or new agents were used.

The existing periodic lower bounds remain -0.64293177476724088 for W=0 and
-0.66063203347329248 for W=1. Their family gaps remain approximately
0.000027128835223944408 and 0.0005974181820038868. Neither numerical limit
is exactly resolved. Refuting these particular mixture extensions does not
refute all mixtures at their energies or change those gaps.

The next concrete step is to integrate the two exact stationary constraints
into the energy and family acceptance paths, propose optimized coefficients,
and require fresh exact PSD and moment replays. Additional overlap conditions
may remain violated afterward. No new coupling-transfer claim is made here;
the previous finite coupling-transfer evidence is preserved. General quantum
representability, generic molecular or long-range/higher-dimensional transfer,
and scalability at requested accuracy remain unproved. The goal stays active.
