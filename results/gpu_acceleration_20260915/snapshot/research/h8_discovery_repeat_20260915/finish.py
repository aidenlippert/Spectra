"""Record this bounded discovery repeat, including its inherited cost boundary."""
from fractions import Fraction as F
import json
import time

from research.h8_discovery_repeat_20260915.run import (
    ROOT, SOURCE, OUT, DISCOVERY, PACKAGE, INPUTS, check_freeze, dump, sha, utc,
)


def run():
    started = time.monotonic()
    protocol = check_freeze()
    read = lambda path: json.loads(path.read_text())
    execution = read(OUT/'execution.json')
    processes = [read(OUT/'runs'/f"{s['name']}.json") for s in protocol['steps']]
    if execution['status'] != 'all_steps_completed' or any(p['status'] != 'passed' for p in processes):
        raise ValueError('Required processes did not all complete; report the failed attempt separately')
    accepted = read(OUT/'exact_replay/complete.json')
    lower = read(OUT/'exact_replay/lower.json')
    if F(accepted['upper_Ha'])-F(accepted['lower_Ha']) != F(accepted['width_Ha']):
        raise ValueError('Exact interval arithmetic mismatch')
    if accepted['target_met'] != (F(accepted['width_Ha']) <= F(1, 625)):
        raise ValueError('Inconsistent target status')
    parent_path = ROOT/protocol['parent_manifest']
    if sha(parent_path) != protocol['parent_manifest_sha256']:
        raise ValueError('Prior milestone manifest changed')
    parent = read(parent_path)
    protected = read(ROOT/'results/h8_spin_import_20260914/preservation_before.json')['files']
    protected.update(parent['files'])
    protected[str(parent_path.relative_to(ROOT))] = {'sha256': protocol['parent_manifest_sha256']}
    parent_seal = ROOT/'results/h8_spin_import_20260914/seal.json'
    if read(parent_seal)['manifest_sha256'] != sha(parent_path):
        raise ValueError('Prior milestone seal mismatch')
    changed = [name for name, data in protected.items() if not (ROOT/name).is_file() or sha(ROOT/name) != data['sha256']]
    if changed:
        raise ValueError(('Inherited work changed', changed))

    singlet = read(OUT/'certificate.json')
    source = read(DISCOVERY/'stage2/export/certificate.json')
    pruned = read(OUT/'runs/prune.log')
    kept = [{**b, 'factor': [r for r in b['factor'] if sum(map(abs, r)) >= 1000]}
            for b in source['core']['blocks']]
    kept = [b for b in kept if b['factor']]
    source['core']['blocks'] = kept
    if source != singlet:
        raise ValueError('Pruning changed content outside the declared row filter')
    old = read(ROOT/'results/h8_spin_import_20260914/local_exact_replay/complete.json')
    if sha(OUT/'certificate.json') == sha(PACKAGE/'certificates/singlet.json'):
        raise ValueError('Expected a newly exported certificate')
    sparse = read(DISCOVERY/'spin_sparse/representation.json')
    normal = read(DISCOVERY/'spin_sparse/build.json')
    stage2 = read(DISCOVERY/'stage2/discovery.json')
    core = {p['name']: p for p in processes}
    total = sum(p['wall_seconds'] for p in processes)
    peak = max(p['peak_child_RSS_bytes'] for p in processes)
    validation = {'frozen_source_and_input_files_verified': len(protocol['source_and_input_hashes']),
                  'previously_protected_files_verified': len(protected), 'changed_previous_files': changed,
                  'original_checker_unchanged': True, 'pruning_matches_declared_rule': True,
                  'new_certificate_differs_from_previous_winner': True}
    audit = {'kind': 'fresh_conditional_H8_discovery_repeated_v1', 'created_UTC': utc(),
             'protocol_sha256': sha(OUT/'protocol.json'), 'trial_count': 1,
             'target_met': accepted['target_met'], 'width_mHa': accepted['width_mHa'],
             'lower_Ha': accepted['lower_Ha'], 'upper_Ha': accepted['upper_Ha'],
             'width_Ha': accepted['width_Ha'], 'prior_milestone_width_mHa': old['width_mHa'],
             'previous_winner_is_preserved': True, 'fresh_optimization_completed': True,
             'inherited_input_names': list(INPUTS), 'cold_from_integrals_completed': False,
             'full_cold_cost_established': False, 'transfer_test_completed': False,
             'competitive_advantage_established': False, 'family_optimum_established': False,
             'factor_rows': accepted['factor_rows'], 'factor_nonzeros': accepted['factor_nonzeros'],
             'singlet_certificate_bytes': (OUT/'certificate.json').stat().st_size,
             'singlet_Gram_entries': sparse['gram_entries'], 'complete_Gram_entries': sparse['total_including_nonsinglet'],
             'coefficient_map_nonzeros': sparse['coefficient_nonzeros'], 'normal_nonzeros': normal['normal_nonzeros'],
             'all_new_process_wall_seconds': total, 'new_discovery_preparation_export_wall_seconds': total-core['exact_replay']['wall_seconds'],
             'exact_replay_process_wall_seconds': core['exact_replay']['wall_seconds'],
             'exact_replay_internal_seconds': accepted['total_replay_seconds'],
             'execution_elapsed_seconds_including_freeze_and_input_checks': execution['elapsed_seconds_including_freeze_and_step_input_checks'],
             'peak_process_RSS_bytes': peak, 'failed_processes': [], 'processes': processes,
             'validation': validation, 'pruning': pruned,
             'post_run_validation_seconds': time.monotonic()-started}
    dump(OUT/'audit.json', audit)
    status = 'succeeded' if accepted['target_met'] else 'missed the target'
    rows = '\n'.join(f"| {p['name'].replace('_', ' ')} | {p['wall_seconds']:.3f} | {p['status']} |" for p in processes)
    text = f'''# H8 discovery repeat: {status}

The frozen two-stage numerical recipe produced a newly exported, exactly checked
interval of **{accepted['width_mHa']:.12f} mHa**. The original
**{old['width_mHa']:.12f} mHa** milestone remains preserved.
This is a new optimization conditional on the declared old inputs, rather than
another replay of the supplied winning factors. One successful repeat does not
establish reliability across random starts or new Hamiltonians.

The exact new lower is `{accepted['lower_Ha']}` Ha. The rational MPS upper is
unchanged. The complete interval covers the entire fixed-eight-electron sector:
the original checker recomputed the upper, expanded the new singlet squares,
checked the separate nonsinglet proof, and subtracted the original-H spin defect
once. No numerical libraries were imported on that accepting path.

## Procedure fixed before execution

The protocol and hashes were frozen at {protocol['frozen_UTC']}. The supplied
optimizer, export code, pruning rule, and original exact checker were unchanged.
All {len(protocol['source_and_input_hashes'])} frozen source/input files still match.
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
{rows}
| **All seven new steps** | **{total:.3f}** | **completed** |

Preparation, optimization, and export account for
**{audit['new_discovery_preparation_export_wall_seconds']:.3f} s**; complete exact
replay accounts for **{core['exact_replay']['wall_seconds']:.3f} s** as a process
({accepted['total_replay_seconds']:.3f} s internally).
Including protocol freezing and repeated input checks, the enclosing run took
**{execution['elapsed_seconds_including_freeze_and_step_input_checks']:.3f} s**.
Peak process RSS was **{peak/1e6:.3f} MB**. Each process has its own parent
resource measurement; platform-native macOS byte units are respected.

These timings begin with prepared maps, a compact checkpoint, the MPS, and a
nonsinglet proof already available. Their earlier construction and unsuccessful
searches are inherited costs, not included or claimed free. Post-run preservation
checks and report generation are outside the enclosing runtime. There is no
complete from-integrals timing or matched speedup claim.

The supplied 90.8186 s replay and the previous 81.0778 s local replay refer to
different executions of the prior certificate. Both were verification timings.

## Search and exported proof

The search still uses {sparse['gram_entries']:,} singlet Gram entries and
{sparse['total_including_nonsinglet']:,} including the inherited nonsinglet component.
It has {sparse['coefficient_nonzeros']:,} stored coefficient-map entries and
{normal['normal_nonzeros']:,} stored normal-matrix entries in this environment.
This run makes no new discovery-compression claim.

The new accepted singlet proof has **{accepted['factor_rows']} factor rows**,
{accepted['factor_nonzeros']:,} nonzero integer factor coefficients, and
{(OUT/'certificate.json').stat().st_size:,} bytes. The pruning rule removed
{pruned['removed_rows']} tiny rows, with a summed square-norm bound of
`{pruned['discarded_square_norm_bound_Ha']}` Ha. The pruned certificate was checked
from scratch. All other certificate content exactly matches the export.
The separate inherited nonsinglet proof still contributes its own factors and cost.

The stage-two numerical dual minimum eigenvalue is
{stage2['best']['minimum_dual_eigenvalue']:.6g}. It supplies no verified ceiling on
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

[PROCEDURE.md]({SOURCE/'PROCEDURE.md'}) records the mathematical acceptance rule,
the exact cost boundary, and the transfer and fragment-consistency diagnostics.
All {len(protected):,} previously protected files checked in this pass are unchanged.

- [Frozen protocol]({OUT/'protocol.json'})
- [Exact complete replay]({OUT/'exact_replay/complete.json'})
- [New accepted singlet certificate]({OUT/'certificate.json'})
- [Measured audit]({OUT/'audit.json'})
'''
    with (SOURCE/'REPORT.md').open('x') as stream:
        stream.write(text)
    files = {}
    for base in (SOURCE, OUT):
        for path in sorted(base.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path not in (OUT/'manifest.json', OUT/'seal.json'):
                files[str(path.relative_to(ROOT))] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
    dump(OUT/'manifest.json', {'kind': 'sealed_H8_conditional_discovery_repeat_v1', 'created_UTC': utc(),
                              'parent_manifest': protocol['parent_manifest'],
                              'parent_sha256': protocol['parent_manifest_sha256'],
                              'files': files, 'file_count': len(files),
                              'milestone': {'fresh_conditional_optimization_target_met': accepted['target_met'],
                                            'cold_from_integrals': False, 'transfer': False}})
    if any(sha(ROOT/name) != row['sha256'] for name, row in files.items()):
        raise ValueError('New evidence changed during sealing')
    seal = {'manifest_sha256': sha(OUT/'manifest.json'), 'new_files_verified': len(files),
            'previously_protected_files_verified': len(protected), 'all_hashes_match': True}
    dump(OUT/'seal.json', seal)
    print(json.dumps({'target_met': accepted['target_met'], 'width_mHa': accepted['width_mHa'],
                      'elapsed_seconds': execution['elapsed_seconds_including_freeze_and_step_input_checks'],
                      'report': str(SOURCE/'REPORT.md'), 'seal': seal}, indent=2))


if __name__ == '__main__':
    run()
