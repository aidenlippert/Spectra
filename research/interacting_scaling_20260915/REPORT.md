# Interacting certificate scaling

The new direct construction independently certifies the original interacting
H8 and H10 Hamiltonians below 1.6 mHa without constructing the complete global
cubic coefficient map. It uses overlapping local spin-adapted blocks and coherent
cross-cluster cubic blocks. This establishes a reusable construction on these
inputs. It has not established the proposed threefold complete-time improvement,
fivefold proof-preparation memory improvement, or general correlation-dependent
scaling.

The frozen rule also passed changed-geometry water and an H4 orbital-space
expansion. The fully coupled H8 path passed. H12 did not reach the target:
its best independently replayed continuation interval is **12.590173 mHa**.
The campaign outcome is a directly constructible, reusable proof family with
demonstrated finite-model scope, accompanied by measured performance limitations.

## Fresh complete controls

| Original input | Direct certified width | Direct fresh time | Direct selected Gram entries | Preserved-family Gram entries |
|---|---:|---:|---:|---:|
| H8, STO-3G | 1.030971387 mHa | 463.317 s | 181,268 | 343,424 |
| H10, STO-3G | 1.557099696 mHa | 1,837.661 s | 702,332 | 1,287,700 |

Both direct runs start with fresh integrals and product-state MPS discovery.
Their clocks include orbital construction, the separate nonsinglet proof, every
rejected proof level, and exact original-H upper/lower replay. The original
Hamiltonian, global fixed particle-number domain, spin allowance and orbital
coefficient allowance are checked. No winning full-family factors or enumerated
state are inputs to these direct runs.

The certified object is the exported rational electronic Hamiltonian. Integral
generation and its coefficient-rounding allowance are recorded, but the original
floating-point integral evaluation is not itself certified. The geometry units
are angstroms; the hydrogen-chain controls use 1.4-angstrom spacing. Nuclear
repulsion shifts both electronic endpoints equally when a total energy is needed.
The frozen geometries, basis choices and source hashes are retained in the
machine-readable protocols. These guarantees do not cover physical-model error.

The H8 reference finished in **324.278 s** with a **0.227658227 mHa** interval.
The direct representation was smaller but slower. Peak proof-preparation RSS
was **446.3 MB direct versus 427.4 MB reference**. On H10 the corresponding
measured peaks were **818.9 versus 789.7 MB** in the initial matched comparison.
Fewer Gram entries therefore did not imply lower peak construction memory.

The first H10 direct level returned a valid **2.467969237 mHa** interval and
missed the target. Its full cost remains in the successful cold clock. Adding
four-orbital local and collective supports, using an exact-checked embedding
of the prior proposal, produced the successful second level. The independent
checker still exported and reconstructed that new proposal from scratch.

The initial H10 full-family reference used 450 seconds of optimization and
returned **10.754481671 mHa**. It stopped just before the common optimizer's
automatic residual-refinement transition. A fresh 840-second main-solve
allocation was consequently declared, under the same 3,600-second whole-pipeline
cap. The failed 450-second attempt remains charged. A finite search that has not
used its complete envelope is not evidence of a reference-family limitation.

The extended fresh reference subsequently passed exact replay at
**0.479085103 mHa in 1,284.551 seconds**, including a main solve that stopped
successfully after 606.102 seconds. Its preparation peak was **815.3 MB**.
Thus the direct H10 method was about **1.43 times slower**, despite using
45.46% fewer selected Gram entries. Both matched sizes currently favor the
optimized reference in complete time. This campaign has not demonstrated a
whole-pipeline performance advantage.

The selected-entry figure describes the accepted representation. The two H10
main-proof preparations together processed **1,117,888 declared Gram entries**,
versus **1,287,700** in the successful reference's one preparation: a **13.19%**
reduction in that cumulative count. The separate magnetic screen and all actual
compute costs remain included in the complete clocks. Summed declared entries
are not peak memory or an arithmetic-work estimate.

The reference is the preserved successful canonical construction with its
established index and spatial masks. It is not the unrestricted cone of every
cubic polynomial. The direct local family is not assumed to be its subset after
orbital rotation. The matched comparisons concern the same original Hamiltonian
and complete accuracy target, rather than identical proof cones. On H8, the
direct family has 34,133 coefficient rows versus 28,461 for the reference,
despite fewer Gram entries; entry count alone does not describe its allocations.

The retained witness bundles were also larger: **5.81 versus 4.62 MB** on H8,
and **15.34 versus 11.53 MB** on H10, direct versus reference. These counts include
the main proof, nonsinglet proof, MPS upper and orbital transform; shared
Hamiltonian descriptions and checker code are excluded. Main exported factor
rows were 286 versus 227 on H8 and 543 versus 322 on H10. The demonstrated
reduction is in selected optimization Gram entries, not exported proof bytes,
complete time or peak memory. [TABLES.md](TABLES.md) reports these quantities
separately.

## What changed mathematically

The direct proof has the form

\[
 H_r-bI=\sum_C O_C^\dagger Q_C O_C+\sum_\alpha B_\alpha^\dagger B_\alpha
       +\mathcal I+R,\qquad Q_C\succeq0.
\]

Global quadratic channels remain present. The collective cubic dictionaries
initially include every pair-supported word, with all cross terms between those
words. Three-, four- and five-orbital windows can enlarge those same coherent
blocks. Local number multipliers multiply the **global** number identity; no
fragment particle number is fixed. All original Hamiltonian terms enter the
coefficient construction or the exact reconstructed remainder allowance.

The maps are generated from those declared words and ideals directly. They are
not slices of a previously built full cubic map. Earlier blocks and ideal terms
remain available at each nested level, so the mathematical optimum cannot get
worse under enlargement. A finite numerical search can still return a worse
candidate; the algorithm retains the best independently accepted interval.

The frozen adaptation rule expands overlapping orbital windows after a failed
complete energy target. It is a geometry-based nested rule, not a learned
selector that has already identified the smallest necessary correlations.
The MPS supplies guide moments and an upper-bound witness; it does not restrict
the unknown physical ground state.

At fixed window width, the representation has O(s^4) collective/quadratic Gram
entries plus O(s k^6) local entries for s spatial orbitals and width k. There is
no proof that fixed k preserves 1.6 mHa as s grows. Indeed, the successful H10
selection uses about 3.9 times H8's entries after enlarging its windows. Those
counts are neither measured runtime exponents nor a scalability theorem.

The complete derivations, residual convention, orbital transformations,
singlet trace and charge-convolution control are in
[MATHEMATICS.md](MATHEMATICS.md). Established related work is distinguished from
this implementation in [SOURCES.md](SOURCES.md).

## Orbital and component controls

The H4 complete-family coordinate control attained **0.000894283 mHa canonical**
and **0.004750749 mHa local**. Preparation improved from 3.657 to 2.995 seconds,
while replay increased from 1.182 to 1.557 seconds. This conditional experiment
does not establish a whole-pipeline orbital speedup. The complete positive/ideal
cone transforms equivalently; the native coefficient-L1 residual penalty need
not. Canonical index restrictions and parity masks are explicitly not copied
into dense local coordinates.

Exact spin-pattern reuse reduced one identical H10 preparation from **320.423
to 71.953 seconds**. All maps and prepared numerical inputs compared
bit-identically. This **4.45-fold preparation-stage improvement** is shared with
the matched reference. It is not a 4.45-fold complete solver improvement.

A direct quadratic M_S=1 proof supplied the H10 nonsinglet screen in **57.85
seconds** for preparation, solve and exact checking. Its lower endpoint was
about 15.45 mHa above the retained ground upper, making it sufficient. The
frozen procedure expands this component when quadratics are insufficient.

See [DEVELOPMENT.md](DEVELOPMENT.md) for conditional attempts and the distinction
between reused discovery inputs and fresh complete calculations.

## Frozen transfer results

Changed-geometry water passed with an original-H interval of
**0.084534412 mHa in 182.389 seconds**, using **111,356 selected Gram entries**.
This is a fresh complete run of the frozen local/collective rule on a non-chain
molecule. It establishes transfer on that declared finite model, without a claim
of experimental accuracy or superiority to existing chemistry methods.

H4 with **6-31g**, expanding four spatial orbitals to eight at the preserved
geometry and particle number, passed at **0.156525949 mHa in 328.437 seconds**,
with **181,268 selected Gram entries**. This is a different Hamiltonian from the
STO-3G fixture. It is an orbital-space expansion test, not another geometry of
the same finite model.

The frozen H12 attempt stopped without a complete interval after **1,502.800
seconds**. Fresh state construction passed in 758.618 seconds and exact upper
checking in 76.601 seconds. All three short magnetic-sector searches missed
the accepting threshold, so the main singlet proof was not attempted. The
quadratic screen's lower lay 10.476684 mHa below the fixed upper; that is not
a full ground-energy interval because the singlet sector was not yet bounded.
The larger screen searches still had substantial reconstruction error.

The predeclared allocation fallback did not trigger: no row/Gram allocation
limit caused this failure. The H12 result is a failure of the frozen complete
procedure under its declared search allocations, not a family-limit theorem.
A separately declared follow-up resumed saved magnetic iterates without changing
that family's operators or exact constraints. Two 300-second solve slices plus
their exact checks took **645.226 additional seconds** and left the screen
**12.716161 mHa** below the retained upper after the spin allowance. This was
still insufficient.

The optimizer had not entered its existing mu=0.03 residual phase. A further
declared diagnostic entered that phase directly, retaining the same maps,
operators and exact constraints. It produced an exactly checked magnetic lower
**0.360475 mHa** below the upper after the spin allowance, with a **0.018093 mHa**
reconstruction allowance,
in **157.226 additional seconds** including initialization and exact checking.
This resolves the magnetic-screen failure within the existing local-cubic family;
it does not by itself bound the singlet ground energy.

The subsequent main-proof initializer exposed a handoff error: the failed frozen
base had no accepted nonsinglet file. That 0.129-second failed attempt is retained.
The corrected driver supplies the newly checked certificate through the existing
explicit source argument and uses a new case directory. It does not modify the
frozen constructor, original failed result or earlier artifacts. The main H12
proof and exact replay remain separately charged within the initial diagnostic
budget. That attempt produced a valid original-model interval of
**452.545612138 mHa** in **1,304.418 seconds**, including **592.309 seconds** of
exact replay at a **1.878 GB** peak. It did not meet the accuracy target.
The initial H12 diagnostics together took **2,107.001 additional seconds**;
the failed frozen run's 1,502.800 seconds remain additional.

One final same-family residual-phase continuation was then declared with a
separate 960-second cap, explicitly beyond that initial diagnostic allowance.
It reuses the same maps and state, changes no operators or accepting equations,
and requires a fresh complete original-model replay. Its result is separate from
the frozen transfer test and is not a matched-resource success claim.

That final continuation returned **12.590172569 mHa**, with exact original-model
replay. It took **906.084 additional seconds**, including **596.964 seconds** of
verification. All H12 diagnostic parent clocks sum to **3,013.085 seconds**, plus
the failed frozen run's **1,502.800 seconds**: **4,515.885 seconds** of recorded
parent work in total. This sum is not a newly timed cold run; intervening work
on other cases and editorial time are excluded. The 1.6 mHa H12 target remains
unmet. The magnetic family was sufficient after refinement, but this main-proof
search does not distinguish an inadequate family from inadequate optimization.

The H12 main preparation uses **1,298,724 Gram entries** but **320,543 coefficient
rows** before spin projection, with 94,225 independent invariant rows. It needs
the separately declared 400,000-row diagnostic envelope; the frozen initial
180,000-row envelope would not accommodate this width-4 preparation. The original
frozen failure occurred earlier, in magnetic screening, so its allocation fallback
had not been triggered. This is another reason to report coefficient-map size
alongside Gram size.

H16 was not allocated because the predeclared H12 success gate was not met.
That resource decision is not evidence of mathematical impossibility at H16.

## Coupled-molecule diagnostic

The H8 path partitions the local orbitals into four two-orbital fragments and
scales every between-fragment term together. At zero coupling the complete
direct interval was **0.004003898 mHa**; at quarter coupling it was
**0.018295612 mHa**, at half coupling **0.059878041 mHa**, and at full coupling
**1.031020918 mHa**. The last interval includes the exact orbital allowance on
the original molecular Hamiltonian. The complete continuation path took
**1,929.606 seconds**; its source integrals and orbital construction are additional.
Between zero and quarter coupling, the retained MPS grew from **48 to 2,449
nonzero integer entries**, and its largest spin-orbital bond from 4 to 64.
Every point used the
same first-level selection of 181,268 Gram entries. This is an
interacting success of a reusable family, not evidence that the constructor
already selects only the smallest necessary couplings.

An exact control freezes the zero-coupling state and evaluates it at full
coupling. Its energy exceeds the original-model ground energy by
**[236.801754036, 237.832777837] mHa**. The allowance includes both orbital
transformation error and the complete ground interval. Thus the coupled state
recovers substantial energy unavailable to this frozen disconnected state.
This comparison includes all variational readjustment; it does not isolate
entanglement energy or prove that a particular SOS block is indispensable.
The control took 1.128 additional seconds and reuses the explicitly charged
source state and ground certificate.

The smaller H4 coupling path also passed all four points, ending at
**0.012779484 mHa** on its original model. Its 308 between-fragment terms include
212 terms that change fragment charges. Independent zero-coupling fragment
controls sum over every allowed local charge and magnetic sector. They use
small local Fock enumerations and charge convolution, and explicitly refuse
nonzero coupling; their success is a control, not the interacting result.

The quarter-coupling solver spent its full 450-second numerical allowance even
though earlier numerical candidates appeared comfortably within the requested
energy width. Its additional residual stopping criterion had not been met.
Those earlier candidates were not independently accepted at that time, so
no hypothetical earlier acceptance is counted as a measured speedup.

## Exact diagnostic and physical-model branch

The H4 singlet-trace repair for the exported pair-supported diagnostic family
imposed 971 exact equations of rank 140, including
27 forced trace-null vectors. A 1/1000 trace mixture passed the independent
exact ideal, Gram positivity and residual-box checks. It took 5.84 seconds
for export, repair and checking, without full particle-sector enumeration.
Its interval floor relative to the frozen upper is only **0.250415 mHa**.
This is a working exact diagnostic, **not an obstruction to reaching 1.6 mHa**.

The CH2 branch retains a defined singlet-triplet energy difference and separates
solver intervals from basis, geometry, correlation and nuclear-motion effects.
The tightly certified minimal model gives roughly 25.49 kcal/mol; higher-basis
numerical controls lie roughly between 10.5 and 13.1 kcal/mol and remain
uncertified physical-model comparisons. No experimental advantage or rigorous
physical error bar is claimed. See [PHYSICAL_MODEL.md](PHYSICAL_MODEL.md).

## Comparisons with other algorithms

Independent exact enumerated controls, freshly constructed from literally
matching rational Hamiltonians, attained approximately **0.800000 mHa** in
**0.861 seconds for H4**, **5.001 seconds for changed-geometry water**, and
**20.559 seconds for expanded-basis H4**. Their clocks include fresh integrals,
independent state/lower discovery and exact verification. They are substantially
faster on these small inputs. Their explicit determinant enumeration is recorded;
it is not a reason to omit the comparison.
The independently constructed exact intervals overlap the corresponding direct
intervals on all three literally matched models.

Separate PySCF numerical FCI controls converged on all five requested models.
Measured process times, including startup, ranged from 0.496 seconds for water
to 7.242 seconds for H12. These runs use the existing pre-rounding spatial
integrals and a tight numerical convergence tolerance. They provide neither
complete fresh-pipeline timings nor rigorous two-sided certificates, and the
convergence tolerance is not a certified energy error. They were run after
the frozen discoveries and did not supply states or proof factors to them.
For the later H12 continuation, the retained checkpoint and rule were declared
separately and do not read the FCI output.

## Reproduction and preservation

[RUNBOOK.md](RUNBOOK.md) gives fresh-run and independent replay commands. All
new artifacts are isolated under this campaign; earlier successful certificates
are preserved. The accepting arithmetic runs independently of the numerical
optimizer. Deliberately enumerated reference calculations are separate controls
and are not inputs to the direct discoveries.

All runs use one heavy local process at a time on the M1 host with 8 GB RAM.
No new cloud instances or external spending were used. Single-run timings do
not estimate run-to-run variance. Source versions, unsuccessful attempts,
complete pipeline clocks and child-stage costs are retained separately.

The complete ledger contains **342 measured child attempts**, totaling
**15,170.247 seconds (4.21 hours)**. This includes unsuccessful searches, tests
and dependency-install attempts, and is not one fresh solver clock. Three process
attempts failed: an initial extension test, the first optional-library build,
and the H12 certificate-handoff initialization. Their corrected follow-ups passed;
the original receipts remain. Accuracy misses are recorded separately from
process failures. The largest measured child RSS was **1.878 GB**. Editorial
work, routine file inspection and final file-integrity accounting are outside
this scientific-stage ledger; parent clocks are not added to child sums again.

Validation includes **27 core regression tests**, **four continuation-input
checks**, exact replay of every claimed interval, an independent 15-state spin
projector comparison, matched-input checks, and overlap with the three independent
exact enumerated controls. The singlet-trace tests cover normalization and the
number/spin ideals through H10, including the H8 singlet dimension 1,764.
Earlier mathematical achievements are retained as references, not overwritten
by this campaign's wider or unsuccessful intervals.
The preservation audit checked **15,892 prior file records** against their
recorded hashes.

[TABLES.md](TABLES.md) contains the measured comparisons and links to individual
receipts. The [accounting record](../../results/interacting_scaling_20260915/accounting.json)
binds all stage receipts by hash. The
[integrity record](../../results/interacting_scaling_20260915/seal.json) binds the
new artifact inventory to the preserved campaign chain. File integrity is
separate from the mathematical checks documented above.
