# Composable response: a second elimination without determinant enumeration

**The first response now admits a second certified elimination on frozen H6.**
Its construction and exact checker use orbital matrices and small local
operator blocks, without generating the retained determinant basis. This
is a reusable mathematical component, not post-hoc compression of the old
retained matrix certificate.

The original **1.000 mHa full interval remains an achieved milestone**.
The frozen fixture is the straight H6 chain at 1.4 Å spacing in canonical
RHF STO-3G orbitals, with the full six-electron sector certified.
All 778 files in its manifest are unchanged, and this pass independently
replayed the complete certificate. The new second response does **not** yet
replace the entire enumerated proof: positivity on the final retained part
remains unproved by the new operator rules.

## Results and what they establish

| Result | Exact evidence | Scope |
|---|---|---|
| Frozen H6 full interval: **1.000 mHa** | Preserved full lower and upper certificate, replayed this pass | Original 1.6 mHa target and compatible old 492-generator family floor already surpassed |
| Frozen H6 second eliminated gap: **0.0402382257 Ha** | Joint occupation-sector bound propagated through the first response, including its residual penalty | New construction and acceptance without full determinant enumeration |
| Second H6 response: order 119, residual penalty **0.000973287 mHa** | Exact Chebyshev scalar inequality | Uniform response error bound; terminal positivity remains separate |
| H6 at 1.6 Å full interval: **1.400 mHa** | New upper and lower independently replayed with the same first-response stopping rule | Complete transfer with explicitly enumerated upper and retained closure |
| H8 second gap estimate: **−0.0083058568 Ha** | Exact lower estimate, so no second positive gap accepted | Bounded larger diagnostic; no new complete H8 interval |
| CH₂ singlet–triplet difference width: **0.0199999583 mHa** | Separate pure S=0/S=1 lower and upper certificates | Exact declared-model result with enumeration; physical accuracy not established |

The retained dimensions are 714 after the first partition and 532 after the
second, computed combinatorially from 924 total N=6 determinants. The new
component constructs neither of those bases. It leaves the operator
obligation K2−eta2 I ≥ 0 on the final 532-dimensional sector. There is no
new compact proof of that last inequality in this pass.

The older fully enumerated reference already had width at most
0.0549943928 mHa. Crossing the old restricted-family floor is a mathematical
milestone relative to that family, not an improvement over every existing
enumerated solver.

## The mathematical change

A separate bound on the second sector's bare energy and on the norm of its
coupling to the first sector was much too loose. Even after improving the
single-sector energy estimates and constructing local CAR coupling bounds,
it produced second-sector lower estimates of −3.599048 Ha on H6 and
−12.365559 Ha on H8.

The successful construction instead certifies the **joint sector** in which
either of the last two spatial orbitals is doubly occupied. Every vector
there has at least two electrons in those orbitals. Density-square tangent
inequalities, a nonnegative occupation multiplier, and a fermionic filling
bound reduce its lower-energy proof to small rational PSD checks.

A joint-sector lower bound carries through the first Schur elimination.
The induced response correction is preserved algebraically; it is not
replaced by a separate coupling-norm penalty when proving the second gap.
Subtracting the first response's certified error still leaves the positive
0.040238 Ha second-sector bound on frozen H6.

The derivation, including both elimination identities and all error shifts,
is in [DERIVATION.md](/Users/aidenlippert/Documents/Spectra/research/composable_response_20260913/DERIVATION.md).
The identity itself belongs to established Feshbach–Schur machinery;
this experiment tests a compact constructible realization.
[Dusson, Sigal and Stamm](https://arxiv.org/abs/2105.02058).

The new orbital PSD matrices are at most 6×6 for H6 and 8×8 for H8.
The accepting path also replays the shared density-tail proof with
21×21 and 36×36 coefficient matrices respectively. Local transition bounds
use seven occupation labels on four spin modes; the other modes remain
operators. The largest local norm majorants have 184 words for H6 and 919
for H8. These are orbital/operator costs, not hidden determinant matrices.

## What the response contributes: matched controls

All controls use the same frozen H6 Hamiltonian, entire N=6 sector, upper,
target b=U−0.001 Ha, first spectral interval and exact factor-verification
machinery. Zero and constant response controls include their actual rounded
residual penalties; coupling is not discarded.

| Construction | Accepted blocks | Failed blocks | Construction and initial check time | Stored proof or obstruction |
|---|---:|---:|---:|---:|
| X=0 | 6/12 | 6 | 0.279 s | 63,109 bytes |
| X=2B*/(M+delta), rationally rounded | 6/12 | 6 | 0.387 s | 111,305 bytes |
| Prescribed polynomial response | 12/12 | 0 | 0.635 s | 654,490 bytes |
| Direct factor of H−b | 12/12 | 0 | 0.443 s | 697,123 bytes |

The shared 924-label block construction takes 0.077 s and is additional to
each row's timer. A separate exact replay verifies all factors and negative
expectations, including complete block coverage. The worst negative
expectations for the tested zero and constant sufficient conditions are
−9.714421 Ha and −5.426109 Ha respectively. These are failures of those
specified sufficient conditions, not molecular ground energies or proofs
against all possible constant responses.

The polynomial response is essential to improvement over these two controls
within the tested closure. When complete determinant matrices are permitted,
a direct factor succeeds and is faster here. Thus the main value of the new
work is its operator-level construction, not an established speed advantage
on H6. The matched mechanisms share dimensions, target and proof machinery;
these timings are not a cost-matched optimization contest over all methods.

## Transfer and larger diagnostic

On H6 at 1.6 Å, the previous HF-based upper improvement was insufficient to
claim a tight full interval. This pass constructed a new upper by diagonalizing
the 200-state conserved block containing HF, after generating all 924 N=6
labels. It then applied the unchanged first-response rule and enumerated
retained-factor machinery at a predeclared width of 1.4 mHa. Both endpoints
passed exact replay. Construction and initial checks took 1.349 s; the saved
closure is 662,080 bytes.

That is a full-interval transfer. It is **not** an enumeration-free transfer.
The new joint-sector estimate on this geometry is −0.0626862054 Ha after
the first penalty, so the second response is not accepted there.

H8 uses the same local partition, tangent construction and stopping rules,
under the same bounded operator budget. It never constructs the 12,870-state
N=8 basis or a many-body matrix. Its second lower estimate is −8.305857 mHa;
therefore no positive second gap or full interval is claimed. H8 had already
been consulted during earlier rule development, so this is a transfer
diagnostic, not a strictly held-out validation of the overall research process.

A negative lower estimate does not establish that the true molecular sector
has no gap. No exact optimality witness for the joint-tangent family was
constructed. Separately, the scalar-envelope approach has exact 2×2
counterexamples consistent with its supplied energy/norm facts. Those prove
that the scalar information alone cannot force positivity; their scope does
not extend to the actual molecular joint sector.

## Costs and inherited work

One cold-process measurement per case and mode gives:

| New joint component | Construction plus initial exact acceptance | Separate exact replay | Discovery peak RSS | Replay peak RSS |
|---|---:|---:|---:|---:|
| Frozen H6 | 4.172 s | 1.264 s | 81.80 MiB | 28.22 MiB |
| H8 diagnostic | 16.099 s | 5.597 s | 141.19 MiB | 76.23 MiB |

These include process startup, imports and frozen input reads. Discovery
includes numerical tangent optimization and its accepting checks. Exact
replay runs with `python -S` and imports no numerical packages. Profiles were
run sequentially with one numerical thread and a 60 s process timeout.
They do not include historical chemistry, density-factor or upper-state
construction. One observation is not a scaling law.

The prescribed new lower component reads Hamiltonian coefficients, the
compact tail certificate and the first response's scalar target/program.
Its final loader does not read any inherited upper amplitudes. Nevertheless,
the complete H6 interval still depends on the frozen upper's 200 FCI-derived
amplitudes, and H8's target was set upstream from 1,000 inherited FCI-derived
amplitudes. Those reference-state discovery costs have not been eliminated.
The fresh geometry's new 200-amplitude upper discovery is explicitly charged.

The input Hamiltonian and orbital choices are ordinary problem data for this
lower-component experiment but were generated by prior chemistry work.
Density factors and certified tails are also inherited, with their prior
cost ledgers linked in the new ledger. No end-to-end fresh-problem runtime
is inferred from the compact component's timings.

There is an additional execution issue. The straightforward composable
recurrence uses 53 base H actions for one first-retained-operator action,
and **12,667 base H actions for one twice-retained-operator action**.
An outer polynomial evaluation on an already supplied right-hand side costs
6,254 H actions. The full K2 count includes constructing its coupling input.
These are formal counts for the current oracle implementation, tested by
exact composition on a small rational example. A molecular K2 action was
not run, and vector-support growth or rational bit growth was not measured.
Small program size therefore does not imply cheap execution.

[The cost ledger](/Users/aidenlippert/Documents/Spectra/results/composable_response_20260913/cost_ledger.json)
records about 132.29 s of measured stages, including failed approaches,
repeated profiles, exact replay and tests. That sum is not total calendar
research time. Ten explicitly enumerated H6 block constructions across the
controls and verification generate 9,240 labels including repeats; they
are kept separate from the zero-enumeration joint-component counts.

## CH₂ pilot

The corrected pilot uses neutral CH₂, one frozen core, and six active
electrons in all six remaining STO-3G orbitals. Its exact S=0 and S=1 proofs
give E_T−E_S approximately between −25.49214 and −25.47959 kcal/mol.
Both upper witnesses are exactly spin-pure, and lower proofs cover the full
declared spin sectors. This requires 400- and 225-determinant matrices.

It validates the model-solving and state-identification path. It does not
establish accuracy against measured CH₂ splitting or an advantage over
conventional FCI, which supplies a faster numerical baseline here. Earlier
incorrect core/electron/spin proposals were invalidated and archived.
The [pilot report](/Users/aidenlippert/Documents/Spectra/research/ch2_validation_20260913/REPORT.md)
defines the geometry, frozen-core treatment, experimental observable,
model limitations and costs.

## Verification and next mathematical obligation

The combined exact replay verifies the preserved 1.000 mHa result, all matched
controls, the 1.400 mHa transfer, all three joint-component outcomes, both
scalar-information obstructions and the CH₂ spin intervals. Fourteen main
focused tests and five CH₂ tests pass, including altered endpoint, PSD,
partition, error-penalty, block-classification and spin-purity refusals.
The final matrix-accounting correction was tested separately; earlier timing
receipts used the new orbital matrix dimension without also naming the shared
tail coefficient matrix. The cost ledger gives the complete dimensions.

Read-only replay from the workspace root:

```sh
cd /Users/aidenlippert/Documents/Spectra
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.composable_response_20260913.verify
```

The original numerical milestone is complete. The next main obligation is
**a compact positivity proof for the final retained operator**, together with
a representation that avoids the straightforward nested action cost. The
joint-sector construction gives a concrete way to carry response effects
through another elimination; it does not yet make those final obligations
free. Further tightening of the already certified original H6 interval is
not the default next objective.
