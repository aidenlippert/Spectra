"""Predeclared H12 resource extension, with the mathematical recipe unchanged."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, dump

PROTOCOL = OUT/'capacity_protocol.json'


def freeze():
    if (OUT/'models/h12_heldout').exists():
        raise ValueError('Declare this fallback before any held-out H12 input is generated')
    h10 = json.loads((OUT/'cold/h10_matched_direct/result.json').read_text())
    costs = [json.loads(p.read_text()) for p in (OUT/'runs').glob('h10_matched_direct_*.json')]
    peak = max(d.get('peak_child_RSS_bytes', 0) for d in costs)
    if not h10['target_met'] or peak > 1200000000:
        raise ValueError('The measured H10 capacity prerequisite is not satisfied')
    dump(PROTOCOL, {'declared_UTC': datetime.now(timezone.utc).isoformat(),
        'base_protocol_sha256': hashlib.sha256((OUT/'heldout_protocol.json').read_bytes()).hexdigest(),
        'adapter_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'trigger': 'An initial H12 preparation explicitly refuses the Gram-entry or coefficient-row envelope',
        'new_Gram_entry_cap': 4000000, 'new_coefficient_row_cap': 400000,
        'H10_measured_peak_RSS_bytes': peak, 'H10_prerequisite_case': 'h10_matched_direct',
        'host_RAM_bytes': 8589934592, 'peak_RAM_is_measured_not_hard_capped': True,
        'fresh_extension_pipeline_budget_seconds': 3600, 'both_attempts_combined_budget_seconds': 7200,
        'same_state_and_proof_algorithms': True, 'same_cluster_levels_and_solve_limits': True,
        'old_attempt_preserved_and_charged': True, 'new_external_spending': False})


def apply(case):
    protocol = json.loads(PROTOCOL.read_text())
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != protocol['adapter_sha256']:
        raise ValueError('Resource adapter changed after its declaration')
    path = case/'design.json'
    design = json.loads(path.read_text())
    if (case/'prepared').exists(): raise ValueError('The envelope may only be set before preparation')
    old = {k: design[k] for k in ('max_Gram_entries', 'max_coefficient_rows')}
    design.update(max_Gram_entries=protocol['new_Gram_entry_cap'], max_coefficient_rows=protocol['new_coefficient_row_cap'])
    path.write_text(json.dumps(design, indent=2)+'\n')
    dump(case/'capacity_extension.json', {'before': old, 'after': {k: design[k] for k in old},
        'protocol_sha256': hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        'operator_supports_or_ideals_changed': False})


def run_conditionally():
    protocol = json.loads(PROTOCOL.read_text())
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != protocol['adapter_sha256']:
        raise ValueError('Resource adapter changed after its declaration')
    if hashlib.sha256((OUT/'heldout_protocol.json').read_bytes()).hexdigest() != protocol['base_protocol_sha256']:
        raise ValueError('The original frozen protocol changed')
    original = json.loads((OUT/'cold/h12_heldout/result.json').read_text())
    messages = ('Declared Gram allocation envelope exceeded before coefficient maps',
        'Declared coefficient row envelope exceeded before map allocation')
    triggers = [str(p.relative_to(ROOT)) for p in (OUT/'runs').glob('h12_heldout*_prepare.log')
        if any(message in p.read_text() for message in messages)]
    if original['target_met'] or not triggers:
        dump(OUT/'capacity_result.json', {'triggered': False, 'original_target_met': original['target_met'],
            'original_seconds': original['wall_seconds']})
        return original
    from research.interacting_scaling_20260915 import cold, pipeline
    frozen = json.loads((OUT/'heldout_protocol.json').read_text())
    for file, value in frozen['source_sha256'].items():
        if hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/file).read_bytes()).hexdigest() != value:
            raise ValueError('The frozen mathematical constructor changed')
    base_execute = pipeline.execute
    def extended_execute(name, steps, deadline=None):
        extended = []
        for label, seconds, command in steps:
            if label == 'prepare' and pipeline.PREFIX+'prepare' in command:
                extended.append(('capacity', 20, [STD, '-B', '-S', '-m', pipeline.PREFIX+'capacity', 'apply', command[-1]]))
            extended.append((label, seconds, command))
        return base_execute(name, extended, deadline)
    old_cold_execute = cold.execute
    pipeline.execute = cold.execute = extended_execute
    try:
        cold.run('h12_capacity_extension', frozen['models']['h12_heldout'], levels=frozen['levels'],
            solve_seconds=frozen['per_level_solve_seconds'], total_seconds=protocol['fresh_extension_pipeline_budget_seconds'])
    finally:
        pipeline.execute, cold.execute = base_execute, old_cold_execute
    result = json.loads((OUT/'cold/h12_capacity_extension/result.json').read_text())
    dump(OUT/'capacity_result.json', {'triggered': True, 'trigger_logs': triggers,
        'initial_attempt_seconds': original['wall_seconds'], 'fresh_extension_seconds': result['wall_seconds'],
        'combined_computation_seconds': original['wall_seconds']+result['wall_seconds'],
        'extension_target_met': result['target_met'], 'extension_is_a_new_fresh_pipeline': True,
        'failed_original_attempt_not_erased': True})
    combined = dict(result)
    combined['wall_seconds'] += original['wall_seconds']
    return combined


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('freeze', 'apply'))
    p.add_argument('case', nargs='?', type=Path); a = p.parse_args()
    freeze() if a.action == 'freeze' else apply(a.case.resolve())
