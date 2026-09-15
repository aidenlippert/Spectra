# H8 accuracy milestone independently verified

**The specified H8 target is achieved.** The supplied certificate independently
replays to **0.767448135463 mHa**, below 1.6 mHa by
0.832551864537 mHa. This records the completed numerical
milestone; it does not replace it with a new requirement.

## Exact result

The new full fixed-N lower is
`-277669368076014353/30000000000000000` Ha, approximately **-9.255645602533813 Ha**.
The unchanged rational MPS upper is approximately **-9.254878154398348 Ha**.
The complete exact fractions and source bindings are in the
[local replay receipt](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/local_exact_replay/complete.json).

| Proof | Complete width, mHa |
|---|---:|
| Previous compact fixed-number construction | 2.849197032649 |
| Preserved stronger large-block comparison | 1.176859795737 |
| New spin-adapted construction, independently checked here | 0.767448135463 |

The new lower improves on the previous compact lower by **2.081748897186 mHa**
and on the preserved stronger lower by **0.409411660274 mHa**.
These comparisons use the same rational Hamiltonian and upper endpoint.
The energies are for the electronic STO-3G H8 model at 1.4 Angstrom spacing;
nuclear repulsion is a separate common constant.

## What was independently accepted

The archive's Hamiltonian, rational MPS, and nonsinglet certificate are byte
identical to the frozen local inputs. All 930 archive-manifest entries match
their supplied hashes. Of 805 packaged source files, 804 match the local files
byte for byte; the additional `portable_h8.py` file was not used by our replay.
The archived source commit matches the local initial commit, `5a0bc6b`.

Our own replay imports the original local checker modules and reads the new
singlet certificate as data. It independently recomputes the actual rational
MPS expectation and expands the new positive squares. It also verifies the
separate nonsinglet proof and charges the original-H spin defect once:

`L = min(L_singlet, L_MS1) - 59/250000000000`.

The result is valid on the **entire eight-electron sector**. Neither NumPy,
SciPy, CVXPY, Quimb, PySCF, Numba, nor SymPy was imported on the accepting path.
The accepting modules and input hashes are recorded in the local receipt.
The new singlet proof needs no residual spectral witness and does not reuse
the stronger comparison certificate to establish its lower bound.

The raw exported b is -9.255633788088 Ha. Its exact averaged coefficient-L1
remainder is `354426294353/30000000000000000` Ha. The checker subtracts this
remainder and the spin defect, obtaining the stated full lower. No floating
optimizer objective, dual optimum, or unverified representation equivalence
is used to accept that energy.

## Why the new representation matters

The successful change retains the full mixed-cubic spin multiplicities. For
each of four charge/parity chains, the magnetic-component dimensions are
112, 368, 368, 112. They decompose into 256 doublet multiplicities and 112
quartet multiplicities: `2*256 + 4*112 = 960`.

The integer raising map U from the +1/2 to +3/2 component satisfies `U U^T=3I`.
The supplied sparse rational kernel basis K has 256 columns. Locally checked
integer identities are `U(6K)=0` and a positive diagonal `(6K)^T(6K)`, whose
diagonal values are 18, 24, and 36. These establish kernel membership and rank
without a numerical rank threshold. All 3,840 mixed words were independently
checked against literal CAR commutators using our local algebra, and the spin
commutator relations passed on all 16 magnetic-component groups.

The old compact point is transported into this larger spin-adapted cone with
a maximum numerical coefficient discrepancy of 8.33e-17.
Two additional random full-Gram checks have relative discrepancies below
3.27e-16. These are numerical
construction diagnostics, distinct from the rational accepting replay.

## Accuracy achieved with a larger search and fewer exported factors

| Quantity | Previous compact attempt | New accepted construction |
|---|---:|---:|
| Singlet Gram matrix entries, sum of squared block sizes | 96,904 | 335,168 |
| Independent symmetric singlet entries | 49,356 | 168,696 |
| Matrix entries including the nonsinglet proof | 115,032 | 353,296 |
| Exported singlet factor rows | 591 | 203 |
| Singlet certificate bytes, including its residual data | 1,901,337 | 1,111,477 |

The new search has 50 blocks, including four 256-dimensional and four
112-dimensional highest-weight blocks. The 203 final factors describe the
accepted proof; they are not the cost or dimension of its discovery problem.
The separate nonsinglet proof adds 316 factor rows and 540,417 bytes.

We verified the pruning exactly: 377 rows were removed and all other
certificate content is unchanged. The sum of their square-norm upper bounds
is `198941/125000000000000000` Ha, about 1.592e-12 Ha. The accepted 203-row
certificate was re-expanded from scratch, independently of that pruning bound.

## Sparse rebuild and numerical portability

The sparse basis and normal system rebuilt successfully from the original
local coefficient maps, without loading the 582 MB dense normal-factor cache.
The imported and locally rebuilt basis arrays are elementwise identical. The
supplied final primal point's selected residual L1 differs locally by only
5.35e-15; the independently measured
normal-map discrepancy is 3.83e-16 relative.

| Numerical storage measurement | Supplied environment | Local reconstruction |
|---|---:|---:|
| Coefficient-map nonzeros | 1,370,610 | 1,494,886 |
| Normal-matrix dimension | 8,533 | 8,533 |
| Stored normal nonzeros | 221,779 | 222,889 |
| Sparse LU factor nonzeros | 954,952 | 954,188 |

The exact storage counts did not reproduce identically. The environments use
different NumPy/SciPy versions, and the local normal matrix contains 31,246
stored entries smaller than 1e-18 in magnitude. Floating cancellation is a
plausible contributor; the exact cause of every count difference was not
isolated. The sparse construction and numerical coefficient agreement were
verified, and the physical interval matches exactly through rational replay.

## Costs and dependencies

The local complete exact replay took **81.078 s**
internally, comprising 32.857 s for the MPS and
48.176 s for the two lower proofs. Its bounded
process took 81.566 s.
All six focused checks passed in 3.530 s. The local sparse
representation/normal build and factorization took 1.367 s
internally, from the existing coefficient maps.

The 3 instrumented local processes total 87.594 process-wall
seconds, with peak recorded process RSS 293.257 MB. Intake, inventory,
light diagnostic probes, and report generation are outside that sum. Local
memory figures come from the platform-aware process runner; the imported
builder's Linux-specific RSS labels are not used for macOS memory accounting.

The external winning trajectory reports a 117.800 s checkpoint and 58.456 s
refinement/export, plus explicitly contracted degree-six MPS moments and other
preparation. We reviewed the winning source path: it uses the old compact
checkpoint, the existing MPS, and full physical coefficient maps; it does not
read the stronger full-cubic factor rows as a teacher. We did not repeat that
optimization or independently total the external exploratory campaign. Its
earlier failed/stopped runs remain preserved in the imported archive.

This establishes the target accuracy and independently checkable representation
and proof. A cold from-integrals runtime, a matched end-to-end speedup, transfer
to a new geometry, and scaling beyond this fixture remain unmeasured here.
No exact obstruction against the old compact spans is claimed.

## Preserved evidence and milestone

All **11,384 previously sealed files** are unchanged. The uploaded
archive and request retain their original hashes, and the imported files remain
unchanged. The archived claim is preserved as supplied; this report records
which parts were independently checked and the local sparsity-count difference.

- [Accepted singlet certificate](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/imported/Spectra_H8_0767448mHa/certificates/singlet.json)
- [Independent complete local replay](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/local_exact_replay/complete.json)
- [Focused representation checks](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/representation_check/receipt.json)
- [Numerical reconstruction comparison](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/representation_check/numeric_bridge.json)
- [Machine-readable audit and achieved milestone](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/audit.json)
- [Supplied derivation](/Users/aidenlippert/Documents/Spectra/results/h8_spin_import_20260914/imported/Spectra_H8_0767448mHa/DERIVATION.md)

**Milestone achieved: a new, independently replayed full H8 interval below
1.6 mHa, with the stronger comparison proof excluded from the accepting path.**
General molecular scalability and optimality of the retained space are separate
unresolved questions. No further optimization campaign was started by this audit.
