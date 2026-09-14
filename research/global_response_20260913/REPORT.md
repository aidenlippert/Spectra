# Spectra: simultaneous response and complete-certificate campaign

This pass produced a substantially cheaper response and exact compact upper
controls. It did **not** discover the inexpensive terminal proof needed for a
new complete solver. The original H6 **1.000 mHa** and transferred H6
**1.400 mHa** milestones remain intact.

The strongest matched execution result is **12,667 → 237 Hamiltonian actions**,
a **53.45-fold reduction**, with the same target and a smaller total response
error allowance. One physical H6 vector took **0.1279 → 0.00644 seconds** on
the same float64 backend, a **19.86-fold observed speedup**. Using an inherited
terminal proof to allocate the error budget further reduces the count to
**175**, while keeping the complete **1.000 mHa** interval. That variant spends
more response error and has an additional proof dependency.

The authoritative results are in
[integrated.json](/Users/aidenlippert/Documents/Spectra/results/global_response_20260913/integrated.json),
[execution.json](/Users/aidenlippert/Documents/Spectra/results/global_response_20260913/execution.json)
and [the derivation](/Users/aidenlippert/Documents/Spectra/research/global_response_20260913/DERIVATION.md).

**An important correction to the earlier baseline.** Existing full cubic SOS
lower proofs already use Hamiltonian coefficients and the particle constraint
without enumerating the full fixed-N determinant space, in both discovery
and acceptance. Their large Gram construction and the separate FCI upper
remain expensive. This campaign does not create the first lower proof that
avoids full determinant enumeration. Its unresolved target is a much cheaper,
small-factor terminal construction that works with the response.

The earlier cubic certificates independently replay to reference-assisted
widths of **0.054994 mHa on H6** and **0.750585 mHa on H8**. The new response
does not improve those stronger existing lower bounds. Instead, they provide
a decisive diagnosis and an explicitly inherited way to close the new
response. Their original discovery uses 217,268 and 1,206,464 Gram entries,
respectively. Their coefficient-space residual dimensions are `(12,66,220)`
and `(16,120,560)`, rather than full fixed-N dimensions 924 and 12,870.

**The terminal diagnosis.** The actual old H6 terminal is `K2−η2I`, with
`K2` already formed from `K1−η1I`. The replayed stronger global lower proves
it positive by at least **0.943512 mHa** at the requested target. The new
uniform-response terminal has a certified margin of **0.944045 mHa**.
There is no negative-terminal obstruction at that H6 target.

A separately charged numerical diagnosis enumerated all 924 fixed-N labels
and the 12 conserved blocks. It found only one terminal eigenvalue below
0.02 Ha; the next low direction across blocks is about 0.07229 Ha. This
motivates an exceptional-direction attack, but those numerical vectors are
not used as a compact complement proof. The positivity statement above comes
from exact SOS replay, not the floating spectra.

**What changed in execution.** The 12,667 count is primarily a product of
the two response degrees. Memoization found no repeated physical Hamiltonian
right-hand sides. The exact identity `F_k=c_k p_(2k)` permits a single
recurrence, but does not remove that product. A simultaneous response against
the combined eliminated sector does remove it. The construction includes
the internal coupling in the combined spectral upper endpoint and bounds
external coupling to the same retained space.

| H6 implementation | Base H actions per vector | Response allowance, µHa | Observed seconds |
|---|---:|---:|---:|
| Original nested | 12,667 | 1.494051 | 0.127917 |
| Exactly equivalent factored recurrence | 12,667 | 1.494051 | 0.187665 |
| Simultaneous, uniform budget | 237 | 0.960339 | 0.006442 |
| Simultaneous, inherited-terminal budget | 175 | 455.772281 | 0.004236 |

The first three rows are the matched uniform-error comparison. The last row
uses the complete-energy margin to choose order 87; its terminal margin is
still **0.489233 mHa**. All rows target the same full H6 lower endpoint.
These timings are individual observations, not statistically established
runtime distributions or end-to-end solver speedups.

Every application here acts on a full **200-entry conserved-block vector**.
The matched benchmark constructs 40,000 H entries for that block. The initial
all-block diagnostic constructs 132,004 H entries. Matrix-free execution
would still carry this vector dependency. Response construction and exact
response checking use orbital matrices and local CAR rules instead; they
do not enumerate these many-body labels. The uniform H6 response descriptor
is 6,212 bytes, excluding its shared Hamiltonian, tail and source code.

**The distinct terminal attacks.**

| Attack | Verified outcome | Limit of the conclusion |
|---|---|---|
| Direct Hamiltonian-ranked two-word factors | H6 lower −41.3853 Ha with 256 atoms; H8 64-atom control is also very weak | Valid lower certificates, no useful terminal closure |
| Frozen local factor procedure | H6 full width 4.008235 Ha; fresh 1.73 Å width 4.502815 Ha after exact spectral residual repair | The tested procedure is inadequate; no family-optimality theorem |
| Low-occupation exceptional space and full complement | H6 complement gap 0.040239 Ha; corrected coupling squared bound 7.042048 Ha² | Scalar Schur penalty about 175.007 Ha is useless; it is not a negative physical witness |
| Collective weighted block domination | Internal coupling squared bound falls 2.983654→2.547938 Ha² on H6, and 7.621997→6.479160 Ha² on H8 | A reusable 14.6–15.0% majorant improvement, not terminal positivity |

The local-family raw H6 objective is −10.176026488 Ha, already several
hartrees below the target. Exact residual repair improves the lower by only
0.025884 mHa. On the fresh geometry it recovers only 0.001762 mHa. Thus
coefficient rounding and residual bounding are not the main cause of these
particular weak outcomes. This is evidence about the bounded solver/export
procedure, not a proof that all factors in the family must fail.

The weighted interference test uses at most 2,187 positive diagonal metrics.
It operates on local occupation blocks, without a configuration-space edge
list. It is not inserted into the frozen response-transfer rule after seeing
the fresh geometry. Its internal coupling map also cannot be reused for the
different exceptional-space partition.

**Transfer and H8.** The rule was frozen before generating H6 at **1.73 Å**:
try the union of the last two double-occupation sectors, otherwise use the
last-double sector; keep a uniform 1 µHa response budget. The saved source
hash still matches. Existing 1.6 Å results are validation, not fresh discovery.

| Case | Selected eliminated sector | Certified sector gap, Ha | H actions | Complete result |
|---|---|---:|---:|---|
| H6, 1.4 Å | Two-sector union | 0.040239 | 237 | 1.000 mHa, inherited cubic terminal proof |
| H6, 1.6 Å | Last double occupation | 0.526513 | 53 | Old 1.400 mHa preserved; new transported closure is 1.410 mHa |
| H8, corrected target | Last double occupation | 0.599426 | 87 | 1.000 mHa, inherited cubic terminal proof; terminal-aware variant needs 67 actions |
| Fresh H6, 1.73 Å | Last double occupation | 0.442816 | 57 | Response passes; new direct full interval is 4.502815 Ha |

The 1.6 Å transported closure deliberately lowers the target by 0.010 mHa to
obtain a verified margin for the different response. It still uses the old
enumerated retained proof. The fresh test includes both energy endpoints;
it does not demonstrate transfer of a tight complete certificate. Sector
gaps in this table are not whole-molecule excitation gaps.

H8 also has a precise target obstruction. Replaying an existing better
2,468-amplitude upper gives an expectation of **−2.152789 mHa** for `H−b_old`.
Therefore the old lower target is physically impossible. The old frozen
1,000-amplitude upper lies **3.152789 mHa** above this better upper, so neither
a 1.0 nor a 1.6 mHa interval with that old upper can work. This says nothing
about a negative direction specifically inside the joint eliminated sector.
Its failed joint-gap proof remains inconclusive. The revised comparison uses
the better upper and a distinct lower target; it is not presented as the
same frozen-upper experiment.

**Upper bounds and the CH2 control.** Exact Wick contraction of rational
disjoint Slater rotations now evaluates an upper using only a one-particle
projector and a proven norm. No determinant expansion is used on that path.
The numerical H optimizer has redundant inactive parameters, records
precision-loss termination, and supplies candidates rather than optima.

The supporting response construction uses
`ψ=(I−αQH)|HF>` and exact numerator/norm moments. It uses five excitation
coefficients on H6 and eight on H8. The current implementation also constructs
`H|HF>` and `Hχ`: **140 distinct labels on H6 and 699 on H8** appear across
these intermediate supports. They remain charged, despite the small recipe.

| Case | Slater upper, Ha | Response-lift upper, Ha | Lift gain over HF, mHa |
|---|---:|---:|---:|
| H6 | −6.1268545064 | −6.1390592879 | 13.285182 |
| H6, 1.6 Å | −5.6126991152 | −5.5611876096 | 18.803450 |
| H8 | −8.9824071637 | −8.9919990234 | 9.591860 |
| H6, 1.73 Å | −5.3693691581 | −5.2440921163 | 23.112517 |

These compact uppers remain about 194–263 mHa above the corresponding tight
references when the better compact candidate is selected. They do not replace
the FCI upper in the tight intervals reported above.

For CH2, common zero-rotation orbitals give exactly pure `S=0` and `S=1`
controls. They are not optimized CH2 results. Replaying the existing model
lower certificates and pairing them with these compact uppers gives widths
**60.358005** and **48.343035 mHa**. With convention `Δ=E_T−E_S`, their gap
interval is **[−0.100972310, +0.007728730] Ha**, which does not resolve the
ordering. The earlier tight CH2 model certificates are preserved. No
experimental comparison or physical-model accuracy is claimed in this pass.

**Cost and verification.** The [cost, dependency and error ledger](/Users/aidenlippert/Documents/Spectra/results/global_response_20260913/accounting.json)
retains measured stages and explicitly marks missing costs. From existing Hamiltonian/tail inputs, uniform
response construction took 2.977 s on H6, 1.907 s on 1.6 Å H6, 6.783 s on H8,
and 2.539 s on fresh H6. Fresh integral/HF, FCI reference and tail generation
together took 1.315 s; this explicitly includes a 400-determinant FCI solve.
The frozen local lower constructions took 13.893 s and 17.574 s.

The final integrated exact run took **57.369 s** with **281,427,968 bytes** peak RSS.
It loaded no NumPy, SciPy, CVXPY or PySCF. It rechecks the large inherited SOS
proofs, both reference uppers, the old 1.6 Å closure, the fresh full interval,
compact uppers and the CH2 controls. Its memory/time therefore include those
dependencies; they are not a cost claim for a small independent terminal
certificate. Historical full-cubic discovery cost and FCI upper discovery
remain separate from this frozen-input replay.

All **27 collected tests pass**, including exact nonunitary-metric refusal,
operator Chebyshev identities, mutated-certificate refusal, local CAR
cross-checks, and independent small determinant-energy oracles for both new
upper paths. The exact integrated run supplies the complete acceptance check.

Some early abandoned branch attempts did not retain reliable elapsed counters.
The ledger marks those costs as unavailable, not zero. Consequently this pass
does not claim an audited total campaign runtime or an end-to-end speedup.
The invalid zero-coupling result, the scaled execution comparison and the
mischaracterized inherited-SOS audit are excluded from all accepted claims.

All **9,486 pre-existing research/results/experiment/document files in the preservation inventory**
match their initial hashes. The prior 855-file manifest is unchanged. This
includes both H6 milestones, the elimination certificates and CH2 proofs.
No publishing, cloud provisioning or purchases were performed.

Replay from the repository root:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -B -S -m research.global_response_20260913.complete
```

The default command is read-only. `--build` is reserved for regenerating this
campaign's derived artifacts. The remaining mathematical target is explicit:
replace the expensive inherited terminal SOS construction while retaining
the new response's measured execution advantage and a useful complete bound.
