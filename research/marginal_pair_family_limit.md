# Certified pair-family limits and the next consistency obstruction

The enlarged pair-transfer family now has independently accepted exact ceilings.
Its optimal periodic lower certificate is bracketed within `1.59e-6` per site
for W=0 and `2.02e-6` for W=1. These are bounds on the fixed certificate family;
they do not establish a physical quantum extension or exact optimum attainment.

The authoritative new directories are
`results/marginal_graded_hubbard8/pair_transfer/<case>/family/`.
`family_summary.json` at the pair-transfer root records exact fractions,
current-source receipt paths and validation evidence.

| Target | Accepted periodic lower/site | New family ceiling/site | Certified gap/site |
|---|---:|---:|---:|
| W=0 | −0.642978716055446 | −0.642977132736826 | 0.000001583318620 |
| W=1 | −0.660850491657565 | −0.660848477151394 | 0.000002014506172 |

Each ceiling uses 81 positive physical source vectors. All 89 retained equations
and the fixed fidelity constraints hold exactly. The lower certificates are
unchanged from the preceding pair-transfer experiment and were freshly replayed
in the new directories. No extra energy improvement is claimed in this turn.

## What the ceiling covers

FAMILY v10 adds the four exact zero pair-transfer moments to every preceding
FAMILY v9 condition. It covers the fixed U=4, t=1, V=1/2, W target, projector
sources, ratio, ceilings and nine sparse shapes. All reflected mean-correct
nearest and range-two profiles, and unrestricted coefficients of the included
charge, hopping, spin, one-spectator and pair-transfer corrections are covered.

The positive source mixture has trace one, satisfies the fidelity inequalities
and annihilates all allowed correction directions. Taking its expectation of
any feasible local lower inequality therefore bounds the lower certificate by
the reported physical local expectation divided by five. That weak-duality
argument gives a ceiling for the certificate family. The mixture need not be
globally extendible, so this is not a physical ground-energy upper bound.

The source cap increases from 85 to 89 only in the explicit pair-family mode.
Old modes retain their old caps, and pair mode requires the full preceding
spectator hierarchy. Fixed pair coefficient fields are still refused in family
certificates: the family verifier checks moments, not a supplied coefficient
choice. The matching driver accepts pair-energy certificates only with the
enlarged cap and still refuses them with older family versions.

Tests verify a stationary mixture, refusal of both old violated mixtures when
relabelled as v10, source-cap boundaries, hierarchy/type checks and the old
driver's refusal path. ENERGY v14 and its physical gates were not changed.

## Construction

The saved last optimizer points were not valid dual certificates. Their Gibbs
mixtures had maximum moment residuals about `0.01465` for W=0 and `0.00636`
for W=1, including positive fidelity residuals. They supplied untrusted vectors
only. Eight dominant atoms were used for W=0 and nine for W=1.

Rounded atoms, deterministic perturbations and all 4096 determinant anchors
produced 4932 and 4980 physical candidates. The numerical LP retained 87 equality
rows and two fidelity inequalities, with zero pricing rounds. It selected 81
source columns in each case. A bounded exact rectangular solve checked all 89
equations and nonnegativity before exporting the proposals. Both passed the
independent production family replay directly; no residual completion, dropped
equation or clipped negative weight was needed.

The new discovery solver admits at most 89 rows. The existing fraction-free
completion helper and its 85-row bound were not changed. The numerical runs
took about 35 seconds each on the local CPU; these observations are not
controlled performance benchmarks or a scalability result.

Fresh independent probes confirm exact zero for all four pair moments, fourteen
one-spectator moments, four spin moments and the unconditioned range-three
hopping moment. The averaged signed-charge law again has an exact stationary
classical Markov extension. These facts do not extend the quantum density matrix.

## Thirty further independent consistency conditions fail

Both new mixtures violate all 30 tested two-spectator charge-hopping conditions.
For distinct sites `i,j,k,l` in a five-site window, consider

\[
C=q_k^p q_l^q B(i,j),\qquad
B(i,j)=\sum_\sigma(c^\dagger_{i\sigma}c_{j\sigma}+c^\dagger_{j\sigma}c_{i\sigma}).
\]

Here `k<l`, and both powers belong to `{1,2}`. Odd-distance hopping uses
`(p,q)=(1,1)` or `(2,2)`; even-distance hopping uses `(1,2)` or `(2,1)`.
These are the PH-even choices. Subtracting the reflected operator gives an odd
five-site `Y`; `T=Y_left−Y_right` is its six-site telescope.

There are 60 parity-allowed primitive products and 30 reflection pairs. This
enumerates the nonconstant-in-both-spectators charge products in this specified
hopping class. It is not an enumeration of all quantum operators.

Independent full-Fock CAR checks establish Hermiticity, reflection, PH and
spin-flip symmetry. Tests reconstruct every action through direct occupation-bit
swaps. Translated coefficients cancel exactly. A two-site row-sum calculation
bounds `||B||` by two; the disjoint charge powers have norm at most one. The
four-term triangle inequality gives `||T||<=8`, without claiming this is sharp.

An exact restricted Gram calculation proves independence modulo the preceding
family. Restrict to off-diagonal two-bit changes with total particle number at
most four. All fixed projectors vanish there; diagonal, spin-exchange and pair
operators have no such entries. The three unconditioned hopping directions
(including two nearest-profile directions) and fourteen one-spectator directions
have rank 17. Adding the new directions raises the exact rank to 47. Thus all
30 are independent modulo the preceding operator class.

| Target | Violated directions | Largest absolute moment | Certified normalized violation using norm bound 8 |
|---|---:|---:|---:|
| W=0 | 30/30 | 0.000547192605073 | 0.0000683990756341 |
| W=1 | 30/30 | 0.00145972476973 | 0.000182465596216 |

The strongest W=0 label is `(i,j,k,l,p,q)=(0,3,2,4,2,2)`; the strongest W=1 label
is `(0,2,3,4,2,1)`. Every reported moment is stored as an exact nonzero rational.
Because the tested operators are symmetry invariant, the violations survive the
stated averaging. They exclude a stationary quantum extension matching these
particular averaged local density matrices.

These violations do not prove that every optimal mixture fails the conditions,
nor that adding the directions must strictly improve the energy. No energy
certificate containing the two-spectator terms has yet been constructed.

## Physical intervals, transfer and remaining scope

The freshly replayed million-site open-chain intervals remain:

| Target | Lower/site | Physical upper/site |
|---|---:|---:|
| W=0 | −0.642981216055446 | −0.610676347051188 |
| W=1 | −0.660854991657565 | −0.618424482369328 |

The full 4096-state local check, overlapping-projector bounds, exact 24-site
physical contraction and 160-bit upper enclosure checks were retained. The
earlier frozen transfer to U=5, V=1/4, W=−1/5 is unchanged; no new coupling or
geometry transfer was attempted here.

The full regression passed 914 tests and 102 subtests in 616.22 seconds. The
focused family/version/refusal checks passed 30 tests. Two subsequently added
two-spectator tests passed separately in 5.47 seconds. All processes are terminal.
Sixteen current receipts and 508 source-hash references were audited. The 114
historical receipts have 3564 unchanged hash references and 149 references that
match preserved pre-change sources; they are not presented as fresh replays.
The exact results and provenance are recorded in `family_summary.json`.

The enlarged-family numerical limit is certified to the finite precision above.
Integrating the next 30 directions remains open. General quantum representability,
generic molecular transfer, long-range or higher-dimensional coverage, and a
requested-accuracy scalability theorem remain unproved. No GPU or paid resource
was used. The broader goal remains active.
