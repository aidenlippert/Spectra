# H8 discovery repeat: succeeded

The frozen two-stage numerical recipe produced a newly exported, exactly checked
interval of **0.725793980368 mHa**. The original
**0.767448135463 mHa** milestone remains preserved.
This is a new optimization conditional on the declared old inputs, rather than
another replay of the supplied winning factors. One successful repeat does not
establish reliability across random starts or new Hamiltonians.

The exact new lower is `-4627801974189358411/500000000000000000` Ha. The rational MPS upper is
unchanged. The complete interval covers the entire fixed-eight-electron sector:
the original checker recomputed the upper, expanded the new singlet squares,
checked the separate nonsinglet proof, and subtracted the original-H spin defect
once. No numerical libraries were imported on that accepting path.

## Procedure fixed before execution

The protocol and hashes were frozen at 2026-09-15T01:58:55.344803+00:00. The supplied
optimizer, export code, pruning rule, and original exact checker were unchanged.
All 693 frozen source/input files still match.
The output directory and Numba cache began empty. MPS moments, their pullback,
the spin representation, and sparse normal system were rebuilt locally.

The first solve used the supplied 130-second budget and penalty 2. The second
used a 180-second budget and penalty 0.03, restarting only from this run's first
stage. It stopped using the supplied numerical criterion. The numerical width
is a selection diagnostic; the accepted width above includes the actual rational
remainder. There was one trajectory and no additional optimizer attempt after
seeing the result. Every instrumented step is listed below.

## Measured cost of this repeat

| Step | Process-wall seconds | Process status |
|---|---:|---|
| mps moments | 4.898 | passed |
| moment pullback | 1.151 | passed |
| sparse build | 1.629 | passed |
| stage1 | 131.964 | passed |
| stage2 | 41.786 | passed |
| prune | 0.076 | passed |
| exact replay | 82.099 | passed |
| **All seven new steps** | **263.604** | **completed** |

Preparation, optimization, and export account for
**181.505 s**; complete exact
replay accounts for **82.099 s** as a process
(81.684 s internally).
Including protocol freezing and repeated input checks, the enclosing run took
**265.367 s**.
Peak process RSS was **677.773 MB**. Each process has its own parent
resource measurement; platform-native macOS byte units are respected.

These timings begin with prepared maps, a compact checkpoint, the MPS, and a
nonsinglet proof already available. Their earlier construction and unsuccessful
searches are inherited costs, not included or claimed free. Post-run preservation
checks and report generation are outside the enclosing runtime. There is no
complete from-integrals timing or matched speedup claim.

The supplied 90.8186 s replay and the previous 81.0778 s local replay refer to
different executions of the prior certificate. Both were verification timings.

## Search and exported proof

The search still uses 335,168 singlet Gram entries and
353,296 including the inherited nonsinglet component.
It has 1,494,886 stored coefficient-map entries and
222,889 stored normal-matrix entries in this environment.
This run makes no new discovery-compression claim.

The new accepted singlet proof has **198 factor rows**,
60,048 nonzero integer factor coefficients, and
1,097,056 bytes. The pruning rule removed
396 tiny rows, with a summed square-norm bound of
`15979/8000000000000000` Ha. The pruned certificate was checked
from scratch. All other certificate content exactly matches the export.
The separate inherited nonsinglet proof still contributes its own factors and cost.

The stage-two numerical dual minimum eigenvalue is
-2.01121e-05. It supplies no verified ceiling on
the best bound expressible by the family. This repeat does not establish why the
old restricted family missed the target or uniquely assign the improvement to
sparsity, added operators, or numerical convergence.

## What changed in the evidence

Fresh optimization from the declared old inputs has now produced another passing
certificate. Previously, local evidence covered only accepting replay of the
supplied winner. Full cold discovery, new-geometry and different-molecule transfer,
larger-size tests, matched comparisons, and physical validation remain open.

The concrete transfer boundary is in the supplied constructor: group indices
42 through 57, fixed H8 paths, a fixed upper, and a residual predictor with
16 orbitals and 8 electrons. A reusable procedure must generate these from model
metadata, rebuild coefficients when spatial symmetry changes, and reconstruct
both upper and lower proofs for the new Hamiltonian. The adaptation rules should
be frozen; block dimensions may grow when justified by accuracy and measured cost.

[PROCEDURE.md](/Users/aidenlippert/Documents/Spectra/research/h8_discovery_repeat_20260915/PROCEDURE.md) records the mathematical acceptance rule,
the exact cost boundary, and the transfer and fragment-consistency diagnostics.
All 12,393 previously protected files checked in this pass are unchanged.

- [Frozen protocol](/Users/aidenlippert/Documents/Spectra/results/h8_discovery_repeat_20260915/protocol.json)
- [Exact complete replay](/Users/aidenlippert/Documents/Spectra/results/h8_discovery_repeat_20260915/exact_replay/complete.json)
- [New accepted singlet certificate](/Users/aidenlippert/Documents/Spectra/results/h8_discovery_repeat_20260915/certificate.json)
- [Measured audit](/Users/aidenlippert/Documents/Spectra/results/h8_discovery_repeat_20260915/audit.json)
