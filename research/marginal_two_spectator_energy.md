# Two-spectator corrections: stronger exact energy certificates

All thirty previously violated two-spectator hopping consistency conditions are
now integrated into the energy verifier and its matching family verifier. Both
new lower bounds strictly exceed the preceding pair-transfer family's certified
ceiling. The W=0 enlarged family is bracketed within 0.00001901 per site; the
W=1 bracket remains much looser. The goal remains active.

The authoritative artifacts are under
`results/marginal_graded_hubbard8/two_spectator/<case>/final/`.
Exact rational values, accepted receipts and source audits are collected in
`two_spectator/summary.json`; `two_spectator/provenance.json` records file hashes.
Decimals below are rounded displays of rational certificates.

| Target U=4, t=1, V=1/2 | W=0 | W=1 |
|---|---:|---:|
| Periodic lower bound per site | -0.6429438983437589 | -0.6607111946158047 |
| Enlarged family ceiling | -0.6429248891913111 | -0.6551630648953223 |
| Certified family gap | 0.000019009152447824 | 0.005548129720482485 |
| Excess over preceding pair-family ceiling | 0.000033234393067114 | 0.000137282535588788 |
| Positive local sources in ceiling certificate | 119 | 118 |
| Million-site open-chain lower bound | -0.6429463983437589 | -0.6607156946158047 |
| Million-site physical variational upper bound | -0.6106763470511881 | -0.6184244823693281 |

The family ceilings bound the best lower certificate in the specified family,
with fixed projector sources, overlap ceilings and sparse shape span. They are
not upper bounds on the physical ground energy. W=0 has a small finite family
gap, without a claim of exact attainment. W=1's numerical limit remains unresolved.
The separate physical energy intervals in the last two rows remain appreciably
wider than either certificate family's optimization interval.

## Physical correction and exact gates

For distinct hopping sites i,j and spectators k,l on five sites, construct
`C = q_k^p q_l^q B(i,j)`, with `B(i,j)` the spin-summed Hermitian hopping.
Reflection-antisymmetrize C to Y, then use the six-site difference
`T = Y_left - Y_right`. Translates telescope to zero on a periodic chain.
Odd hopping distances use powers (1,1) and (2,2); even distances use (1,2)
and (2,1). The sixty primitive products give thirty reflection pairs.
The previous independent CAR probe established symmetry, cancellation and
independence modulo the previous family. These are stationary overlap
consistency constraints; they do not change the target Hamiltonian.

`marginal_two_spectator_hopping.py` implements the thirty canonical actions.
ENERGY v15 adds their exact projected matrix before the local PSD test.
FAMILY v11 requires all thirty expectations to vanish exactly and requires the
complete preceding pair-transfer hierarchy. Fixed coefficients are forbidden
in a family certificate: the family constrains moments instead. Older energy
versions refuse the new field, and the matching driver refuses a v15 energy
with a preceding family cap. The old modes and their source limits remain.

The new mode increases the nine-shape source cap from 89 to 119. Local energy
verification still covers 94 blocks and all 4096 Fock states, with maximum local
PSD dimension 200. Rational weight strings remain within the existing
4096-character bound. The selected longest weights have 3404 and 3710 characters.
There is no enlarged physical window or new opening interaction cost.

Independent CAR reconstruction confirms all thirty moments are exactly zero in
both selected mixtures. The four pair-transfer, fourteen one-spectator, four spin
and coherent hopping moments also vanish. Their averaged signed-charge laws
admit exact classical Markov extensions. None of these checks establishes an
extension of the full quantum density matrix or general representability.

## Discovery and unsuccessful searches

The thermal search uses 118 coefficients, full spectral gradients and curvature
preconditioning, two temperatures and a maximum of 500 full spectral evaluations
per run. W=0 used 469 evaluations; W=1 used 500, followed by a separate 500-evaluation
polish. These were iteration- or budget-limited searches. The exported Gibbs
atoms and floating-point solver outputs are proposals, with no acceptance power.
W=0's initial lower and W=1's polished lower were selected only after exact replay.

Initial sampled family fits returned the trivial ceiling zero. Global column
pricing over all 94 blocks, with up to two eigenvectors per block and forty
rounds, produced 9276 W=0 candidates and 10506 W=1 candidates. Exact selected
column reconstruction gave the accepted W=0 ceiling and an intermediate W=1
ceiling of -0.653901827567574. A separate W=1 run using a radius-0.02 clipped
pricing vector executed forty rounds and remained at zero. Sampling near the
polished lower also returned zero. Those unsuccessful proposals are retained.

Resuming the existing W=1 candidate ledger for forty further global rounds
produced 14630 candidates and the selected ceiling -0.6551630648953223. The last
reduced eigenvalues remained negative (about -1.45 for W=0 and -2.61 for resumed
W=1), so pricing convergence is not established. Numerical LP optimality on a
finite candidate set does not certify the full family optimum. Exact feasible
mixtures do certify the displayed ceilings regardless of that search failure.

The resume helper accepts at most 13000 input candidates and emits at most
20000. The latest 14630-candidate output cannot be fed straight back into it.
A further attempt needs a justified candidate reduction or a revised pricing
strategy; simply repeating the present run is not an established solution.
Observed numerical run times of roughly 33–231 seconds are local CPU timings,
not controlled scalability or GPU benchmarks. No GPU or paid resources ran.

## Frozen transfer and ablation

The W=0 recipe was transferred to U=5, t=1, V=1/4, W=-1/5 on the same half-filled
chain. Onsite and density profiles were scaled to the target and the range-two
profile shifted as required by the new coupling. Projectors, overlap ceilings,
penalties and correction coefficients were frozen. Only the local spectral
threshold was recomputed and then certified exactly.

| Million-site open-chain result | Energy per site |
|---|---:|
| Frozen new recipe lower | -0.5226432283437589 |
| Same recipe with only thirty new corrections removed | -0.5230567883437589 |
| Previous frozen pair-transfer recipe lower | -0.5230701060554458 |
| Physical variational upper | -0.4885616802989547 |

Removing only the thirty corrections loses exactly `10339/25000000 = 0.00041356`
per site. The new recipe improves on the previous frozen recipe by approximately
0.00042687771168696. All three energy recipes passed fresh exact replay and the
comparison checks their frozen fields. These are comparisons between fixed
recipes, without a ceiling for a reoptimized older family at the new target.
This is useful finite coupling transfer; geometry, filling and short range are
unchanged. Generic molecular, long-range and higher-dimensional transfer remain
unproved.

## Validation and remaining work

Focused production and boundary tests passed (44 tests plus one matching test).
The final full suite passed **945 tests and 102 subtests in 849.76 seconds**.
One existing calibration test returns a value and triggers a pytest warning.
The suite includes independent per-label CAR comparisons, full-image projection,
negative-PSD refusal, malformed input and older-version refusals, source-cap
boundaries and rejection of relabeled older violating mixtures. No production
source changed after the full suite was collected.

All selected energy, family, strict separation, independent closure and transfer
receipts were audited against current source hashes. Earlier receipts remain
historical evidence with changed sources checked against preserved snapshots;
they are not represented as fresh replays under the new verifier. Full
construction diagnostics, unsuccessful proposals, candidate ledgers and logs
are retained. All launched computation and validation jobs have terminated.

The next unresolved computational issue is the W=1 family gap. Further quantum
consistency constraints have not yet been tested on these new mixtures. General
representability, exact family attainment and cost at requested accuracy remain
open; this finite-chain progress does not settle general quantum chemistry.
