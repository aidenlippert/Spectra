"""Record the achieved H8 milestone and preserve its independent receipts."""
from fractions import Fraction as F
import datetime
import hashlib
import json
from pathlib import Path

from research.h8_spin_import_20260914.budget import ROOT, OUT

SOURCE = ROOT/'research/h8_spin_import_20260914'
PACKAGE = OUT/'imported/Spectra_H8_0767448mHa'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def run():
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed validation pass')
    read = lambda path: json.loads(path.read_text())
    intake = read(OUT/'intake.json')
    verified = read(OUT/'local_exact_replay/complete.json')
    lower = read(OUT/'local_exact_replay/lower.json')
    tests = read(OUT/'representation_check/receipt.json')
    sparse = read(OUT/'representation_check/sparse_rebuild_receipt.json')
    bridge = read(OUT/'representation_check/numeric_bridge.json')
    supplied = read(PACKAGE/'RESULT.json')
    processes = [read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    if any(p['status'] != 'passed' for p in processes):
        raise ValueError('Incomplete or failed verification process')
    if not verified['target_met'] or not tests['successful']:
        raise ValueError('Required acceptance gate did not pass')
    if F(verified['upper_Ha'])-F(verified['lower_Ha']) != F(verified['width_Ha']):
        raise ValueError('Exact interval arithmetic mismatch')
    if F(verified['lower_Ha']) != F(supplied['lower_Ha']):
        raise ValueError('Imported claim differs from independent acceptance')
    if sha(intake['archive']) != intake['archive_sha256'] or sha(intake['request']) != intake['request_sha256']:
        raise ValueError('Supplied archive or request changed')
    package_manifest = read(PACKAGE/'SHA256.json')
    if any(sha(PACKAGE/name) != digest for name, digest in package_manifest.items()):
        raise ValueError('Imported package modified')
    protected = read(OUT/'preservation_before.json')['files']
    changed = [name for name, data in protected.items() if not (ROOT/name).is_file() or sha(ROOT/name) != data['sha256']]
    if changed:
        raise ValueError(('Inherited work changed', changed))
    L, U = F(verified['lower_Ha']), F(verified['upper_Ha'])
    gain = float(F(verified['lower_improvement_vs_previous_compact_Ha'])*1000)
    strong_gain = float(F(verified['lower_improvement_vs_previous_strong_Ha'])*1000)
    wall = sum(p['wall_seconds'] for p in processes)
    peak = max(p['peak_child_RSS_bytes'] for p in processes)
    audit = {'kind': 'H8_accuracy_milestone_independently_achieved',
             'claim_status': 'Verified by unchanged original local accepting code',
             'lower_Ha': str(L), 'upper_Ha': str(U), 'width_Ha': str(U-L),
             'width_mHa': float((U-L)*1000), 'target_mHa': 1.6,
             'target_met': True, 'target_margin_mHa': float((F(1, 625)-(U-L))*1000),
             'new_lower_improvement_mHa': gain, 'improvement_over_old_strong_mHa': strong_gain,
             'local_exact_replay': 'results/h8_spin_import_20260914/local_exact_replay/complete.json',
             'local_tests_passed': tests['tests_run'], 'original_checker_modified': False,
             'original_input_bytes_identical': intake['input_comparison'],
             'source_comparison_counts': intake['source_comparison_counts'],
             'inherited_files_verified': len(protected), 'changed_inherited_files': changed,
             'imported_manifest_hashes_verified': len(package_manifest),
             'singlet_Gram_matrix_entries': tests['matrix_entries'],
             'singlet_independent_symmetric_entries': tests['independent_symmetric_entries'],
             'complete_Gram_matrix_entries': tests['matrix_entries']+18128,
             'accepted_singlet_rows': lower['singlet']['factor_rows'],
             'accepted_singlet_nonzero_factor_coefficients': lower['singlet']['factor_nonzeros'],
             'sparse_rebuild': sparse, 'numeric_bridge': bridge,
             'local_instrumented_processes': len(processes), 'local_recorded_process_wall_seconds': wall,
             'local_peak_process_RSS_bytes': peak,
             'fresh_numerical_discovery_repeated_locally': False,
             'cold_end_to_end_cost_verified': False, 'general_scalability_demonstrated': False,
             'old_compact_span_impossibility_proved': False,
             'no_next_research_campaign_automatically_started': True}
    write_json(OUT/'audit.json', audit)
    text = f'''# H8 accuracy milestone independently verified

**The specified H8 target is achieved.** The supplied certificate independently
replays to **{verified['width_mHa']:.12f} mHa**, below 1.6 mHa by
{audit['target_margin_mHa']:.12f} mHa. This records the completed numerical
milestone; it does not replace it with a new requirement.

## Exact result

The new full fixed-N lower is
`{L}` Ha, approximately **{float(L):.15f} Ha**.
The unchanged rational MPS upper is approximately **{float(U):.15f} Ha**.
The complete exact fractions and source bindings are in the
[local replay receipt]({OUT/'local_exact_replay/complete.json'}).

| Proof | Complete width, mHa |
|---|---:|
| Previous compact fixed-number construction | 2.849197032649 |
| Preserved stronger large-block comparison | 1.176859795737 |
| New spin-adapted construction, independently checked here | {verified['width_mHa']:.12f} |

The new lower improves on the previous compact lower by **{gain:.12f} mHa**
and on the preserved stronger lower by **{strong_gain:.12f} mHa**.
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
remainder is `{lower['singlet']['residual_l1']}` Ha. The checker subtracts this
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
a maximum numerical coefficient discrepancy of {tests['warm_start_max_coefficient_error']:.3g}.
Two additional random full-Gram checks have relative discrepancies below
{max(tests['random_full_Gram_transport_relative_errors']):.3g}. These are numerical
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
{bridge['absolute_residual_l1_difference']:.3g}; the independently measured
normal-map discrepancy is {bridge['normal_map_relative_error']:.3g} relative.

| Numerical storage measurement | Supplied environment | Local reconstruction |
|---|---:|---:|
| Coefficient-map nonzeros | 1,370,610 | {sparse['coefficient_map_nonzeros']:,} |
| Normal-matrix dimension | 8,533 | 8,533 |
| Stored normal nonzeros | 221,779 | {sparse['normal_stored_nonzeros']:,} |
| Sparse LU factor nonzeros | 954,952 | {sparse['LU_factor_nonzeros']:,} |

The exact storage counts did not reproduce identically. The environments use
different NumPy/SciPy versions, and the local normal matrix contains 31,246
stored entries smaller than 1e-18 in magnitude. Floating cancellation is a
plausible contributor; the exact cause of every count difference was not
isolated. The sparse construction and numerical coefficient agreement were
verified, and the physical interval matches exactly through rational replay.

## Costs and dependencies

The local complete exact replay took **{verified['total_replay_seconds']:.3f} s**
internally, comprising {verified['upper_replay_seconds']:.3f} s for the MPS and
{verified['lower_replay_seconds']:.3f} s for the two lower proofs. Its bounded
process took {next(p['wall_seconds'] for p in processes if p['name']=='local_exact_replay'):.3f} s.
All six focused checks passed in {tests['seconds']:.3f} s. The local sparse
representation/normal build and factorization took {sparse['build_and_factor_seconds']:.3f} s
internally, from the existing coefficient maps.

The {len(processes)} instrumented local processes total {wall:.3f} process-wall
seconds, with peak recorded process RSS {peak/1e6:.3f} MB. Intake, inventory,
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

All **{len(protected):,} previously sealed files** are unchanged. The uploaded
archive and request retain their original hashes, and the imported files remain
unchanged. The archived claim is preserved as supplied; this report records
which parts were independently checked and the local sparsity-count difference.

- [Accepted singlet certificate]({PACKAGE/'certificates/singlet.json'})
- [Independent complete local replay]({OUT/'local_exact_replay/complete.json'})
- [Focused representation checks]({OUT/'representation_check/receipt.json'})
- [Numerical reconstruction comparison]({OUT/'representation_check/numeric_bridge.json'})
- [Machine-readable audit and achieved milestone]({OUT/'audit.json'})
- [Supplied derivation]({PACKAGE/'DERIVATION.md'})

**Milestone achieved: a new, independently replayed full H8 interval below
1.6 mHa, with the stronger comparison proof excluded from the accepting path.**
General molecular scalability and optimality of the retained space are separate
unresolved questions. No further optimization campaign was started by this audit.
'''
    with (SOURCE/'REPORT.md').open('x') as stream:
        stream.write(text)
    files = {}
    for root in (SOURCE, OUT):
        for path in sorted(root.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path not in (OUT/'manifest.json', OUT/'seal.json'):
                files[str(path.relative_to(ROOT))] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
    manifest = {'kind': 'sealed_H8_spin_import_validation_v1',
                'created_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'parent_manifest': 'results/sector_quotient_20260914/manifest.json',
                'parent_sha256': intake['parent_manifest_sha256'],
                'source_archive_sha256': intake['archive_sha256'], 'files': files,
                'file_count': len(files), 'total_bytes': sum(row['bytes'] for row in files.values()),
                'milestone': {'H8_full_interval_below_1p6mHa_independently_verified': True,
                              'general_scalability': False, 'cold_discovery_reproduced': False}}
    write_json(OUT/'manifest.json', manifest)
    if any(sha(ROOT/name) != value['sha256'] for name, value in files.items()):
        raise ValueError('New evidence changed during sealing')
    seal = {'manifest_sha256': sha(OUT/'manifest.json'), 'new_files_verified': len(files),
            'inherited_files_verified': len(protected), 'all_hashes_match': True}
    write_json(OUT/'seal.json', seal)
    print(json.dumps({'width_mHa': verified['width_mHa'], 'milestone_achieved': True,
                      'report': str(SOURCE/'REPORT.md'), 'seal': seal}, indent=2))


if __name__ == '__main__':
    run()
