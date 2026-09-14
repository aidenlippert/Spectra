# Spectator hopping: strict escape and the remaining pair-transfer obstruction

The complete nonconstant one-spectator charge-hopping class gives exact lower
certificates that strictly exceed the previous fixed spin-family ceiling on
all four targets. This rules out obtaining the new bounds by retuning only the
previous family's coefficients. The enlarged family's numerical optimum is
still unresolved: all eight primal search passes hit their iteration limit,
and its certified lower-to-ceiling gaps remain loose.

The authoritative artifacts are in
`results/marginal_graded_hubbard8/spectator_hopping/<case>/polished/`.
`combined_summary.json` contains exact fractions, receipt paths and audit results;
`numerical_provenance.json` identifies each selected dual proposal's seed and
preserved numerical driver. Discovery proposals are nonaccepting until replayed.

## Certified physical intervals and fixed-family separation

The target is the open, half-filled, million-site chain

\[
H=4\sum_i n_{i\uparrow}n_{i\downarrow}
-\sum_{i,\sigma}(c^\dagger_{i\sigma}c_{i+1,\sigma}+\mathrm{h.c.})
+\tfrac12\sum_i q_iq_{i+1}+W\sum_i q_iq_{i+2},
\qquad q_i=n_{i\uparrow}+n_{i\downarrow}-1.
\]

Decimals below are rounded displays of exact rational receipts.

| W | Open lower/site | Physical upper/site | Periodic excess over old spin-family ceiling | Current periodic family gap |
|---|---:|---:|---:|---:|
| 0 | −0.643064190986 | −0.610676347051 | 0.0000623508253 | 0.000209786662 |
| +0.1 | −0.643712090957 | −0.611451160583 | 0.0000568915972 | 0.000220135388 |
| −0.1 | −0.642656056991 | −0.609901533519 | 0.0000985000890 | 0.000146692895 |
| +1 | −0.661010393817 | −0.618424482369 | 0.0000858084851 | 0.000512088956 |

Exact open lower fractions, in the same order:

- `-16076604774649909/25000000000000000`
- `-16092802273917671/25000000000000000`
- `-16066401424763747/25000000000000000`
- `-413131496135603/625000000000000`

Every physical upper was freshly evaluated using the filtered-state recurrence,
including exact 24-site contraction checked against its 160-bit enclosure. The
uppers are unchanged from the preceding milestone.

Each `family_escape.json` freshly replays the prior spin-family witness using
current production code and matches the physical target, projector sources,
ratio, ceilings and nine-shape sparse span. Its strictly positive exact excess
covers unrestricted coefficients in that previous family and its reflected,
mean-correct profiles. A family ceiling limits attainable lower certificates;
it is not a physical energy upper. This comparison does not cover arbitrary
projector sources, larger supports or all certificate constructions.

## Complete one-spectator class

Let \(B(i,j)=\sum_\sigma(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})\).
For distinct spectator \(k\), set \(C(i,j,k)=q_k^p B(i,j)\), with \(p=2\)
for odd hopping distance and \(p=1\) for even distance. On five sites define
\(Y=C(i,j,k)-C(4-j,4-i,4-k)\); on six sites use
\(T=Y_{\rm left}-Y_{\rm right}\). These corrections cancel under translation,
so they add neither a physical interaction nor an opening penalty.

There are 14 canonical reflected classes: nine quadratic-charge directions and
five signed-charge directions. Because a function on \(q\in\{-1,0,1\}\) has
the form \(a+bq+cq^2\), these exhaust the nonconstant single-spectator charge
functions multiplying real hopping after the stated particle-hole and reflection
projection. The three unconditioned odd-distance hopping directions account for
the constant part. This completeness statement is restricted to this operator
class.

`discovery/spectator_hopping_basis.py` independently reconstructs CAR actions
over all 4096 states and checks Hermiticity, spin-flip, particle-hole and
reflection symmetries, and periodic cancellation at sizes 8, 9 and 10. Its
integer Gram matrix on the at-most-two-electron sector has rank 17 for the
three old unconditioned directions plus the 14 new directions. Prior diagonal
and spin corrections have no single-hop off-diagonal entries, while the actual
fixed projector sources have particle numbers 5, 6 or 7 and vanish in this
sector. Thus all 14 directions are independent modulo the prior fixed-source,
mean-correct variable span.

The parent `basis_probe.json` files refer to the old spin witnesses: every one
of their 14 moments is exactly nonzero. The final `one_spectator_overlap.json`
files refer to the new witnesses: all 14 moments are exactly zero. Independent
spin and unconditioned hopping closures also pass. The exact ablation witnesses
show that removing only the spectator term makes the submitted local residual
negative, while restoring it makes that same vector's residual nonnegative.
This ablation proves necessity at those fixed coefficients; the separate family
escape establishes separation from all tuning inside the older family.

## Production and bounded numerical work

ENERGY v13 adds an exact, bounded `spectator_hopping` coefficient dictionary;
older versions reject that field. FAMILY v9 adds all 14 exact zero-moment
constraints and raises the source cap from 71 to 85. Every selected mixture
uses 85 integer-vector physical sources. The family input refuses fixed
spectator coefficients, and the matching replay requires the new family mode.
All previous constraints and refusal paths remain active.

The lower replay retains 94 local blocks spanning 4096 states, maximum local
PSD dimension 200, 151 total PSD checks and the sparse-entry cap of 64. The new
helper builds exact physical actions and projects them with the existing
congruence map. No block-size or asymptotic scalability improvement is claimed.

The primal chart has 85 variables: two hopping profiles, two penalties, nine
sparse corrections, 52 signed-charge corrections, one unconditioned hopping
telescope, four spin telescopes, 14 spectator telescopes and the local lower
value. Both first and polishing passes were bounded to 240 optimizer iterations
and 250 matrix evaluations. All eight passes reached the iteration limit;
actual matrix evaluations were 240 or 241. Fresh physical affine reconstruction
differed from its numerical matrix by at most about `3.56e-15`. Exact replay,
not that floating-point agreement, accepts the lower certificates.

The dual proposal has 85 exact equations and retains all of them in the
rectangular rational solve. Determinant anchors include the nonzero spin
diagonals. Candidate caps remain 1200 sampled states, 4096 anchors and 13000
total after pricing; LP tolerance is `1e-10`, support threshold `1e-14`.

| Target | Selected seed | Eigenvectors per priced block | Pricing rounds | Candidates |
|---|---|---:|---:|---:|
| W=0 | First accepted lower | 1 | 60 | 8130 |
| W=+0.1 | First accepted lower | 1 | 60 | 8008 |
| W=−0.1 | Polished lower | 2 | 40 | 8482 |
| W=+1 | Polished lower | 2 | 40 | 8470 |

The first W=−0.1 dual attempt produced an inconsistent 84-column exact basis
and wrote no proposal. The bounded rerun produced an accepted 85-source witness.
The W=1 rerun improved its ceiling. The seed and number of rounds changed too,
so these observations do not isolate an effect of the eigenvector count.
The selected W=0 and W=+0.1 parent proposals were copied into their polished
directories and independently matched and replayed there. Their preserved
pre-multieigenvector driver is identified in the provenance manifest.

Four initial ablation invocations started before their required energy receipt
existed and exited with `FileNotFoundError`. Those logs are retained as
`ablation_before_energy_ready.log`. After the energy replays completed, all four
ablations were rerun successfully. This was an execution-order error, not a
failed positivity test. No computation remains running at this checkpoint.

## Next exact obstruction: pair transfer

Define \(d_i^\dagger=c^\dagger_{i\uparrow}c^\dagger_{i\downarrow}\) and
\(J(i,j)=d_i^\dagger d_j+d_j^\dagger d_i\). The four canonical pairs
`(0,1)`, `(0,2)`, `(0,3)` and `(1,2)` define reflected five-site differences
and six-site translated telescopes as above. Every one of the four moments is
nonzero in every selected new witness.

| Target | Strongest pair | Moment | Absolute moment / norm bound 4 |
|---|---|---:|---:|
| W=0 | (0,1) | +0.001213418315 | 0.000303354579 |
| W=+0.1 | (0,1) | +0.000886622458 | 0.000221655614 |
| W=−0.1 | (0,1) | +0.001564405034 | 0.000391101259 |
| W=+1 | (0,3) | −0.001289508262 | 0.000322377066 |

`discovery/pair_transfer_probe.py` checks the full-Fock CAR actions,
Hermiticity, particle-number conservation, particle-hole/reflection/spin-flip
symmetries and periodic cancellation at sizes 8, 9 and 10. An independent
two-site calculation gives \(\|J\|=1\), hence the telescope has norm at most 4
by the triangle inequality. This is an upper bound, not an exact norm claim
for each telescope.

On the two-electron basis with one doubly occupied site, the exact Gram matrix is

```
[[ 8, 0,  0, -4],
 [ 0, 8,  0,  0],
 [ 0, 0, 12,  0],
 [-4, 0,  0, 12]]
```

It has rank four. Diagonal corrections have zero off-diagonal entries there;
spin exchange cannot move a double occupancy; single hopping changes two bits
where pair transfer changes four; the fixed projectors again vanish. These are
four independent new directions outside the current family.

The moments persist under particle-hole/reflection averaging. They therefore
exclude a stationary quantum extension matching that averaged six-site density
matrix. They do not rule out every inhomogeneous extension of a single marginal.
No pair-transfer correction has yet been integrated or optimized.

The averaged signed-charge projection nevertheless has an exact order-five
stationary classical Markov extension: all 243 prefix/suffix equalities,
stochastic rows, shift compatibility, stationarity and 729 reconstructed
probabilities pass. Classical charge consistency and the tested coherent
conditions are insufficient for quantum representability.

## Validation and remaining scope

The focused spectator-plus-spin run passed 28 tests in 111.41 seconds. The full
current regression passed 847 tests and 102 subtests in 661.31 seconds. JUnit
records 847 testcase elements and 949 tests including subtests, with zero
errors, failures or skips. One existing warning concerns a test returning a
tuple. Four basis probes and 36 final receipts are accepted; all 1348 recorded
source-hash entries match their current files. The checkpoint collector only
audits existing receipts; it is not a replacement for their mathematical replay.

The current optimization gap remains open, especially at W=1. The four
pair-transfer directions provide a verified next compact extension, but closing
them would still not prove general representability. Transfer to arbitrary
molecular Hamiltonians, long-range or higher-dimensional models, and scalability
at a requested accuracy remain unproved. No GPU was used in this milestone.
The broader goal remains active.
