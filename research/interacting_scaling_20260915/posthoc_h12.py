"""Bounded diagnosis of the H12 magnetic-search allocation, outside the freeze."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.pipeline import execute, PREFIX

PROTOCOL = OUT/'posthoc_h12_protocol.json'


def declare():
    result = json.loads((OUT/'cold/h12_heldout/result.json').read_text())
    if result['target_met']: raise ValueError('This diagnostic concerns the observed failed frozen run')
    dump(PROTOCOL, {'declared_UTC': datetime.now(timezone.utc).isoformat(),
        'frozen_result_sha256': hashlib.sha256((OUT/'cold/h12_heldout/result.json').read_bytes()).hexdigest(),
        'source_screen': 'h12_heldout_base_magnetic_level1',
        'reason': 'The enlarged magnetic searches restarted from zero and ended with substantial reconstruction error',
        'post_hoc_not_a_frozen_heldout_success': True,
        'maximum_additional_seconds': 2400, 'same_family_continuation_slices_seconds': [300, 300],
        'main_proof_if_screen_succeeds': {'widths': [2, 3, 4], 'collective_widths': [3, 4], 'solve_seconds': 450,
            'Gram_entry_cap': 4000000, 'coefficient_row_cap': 400000},
        'all_prior_discovery_and_failed_attempts_are_additional': True,
        'source_state_and_maps_reused_explicitly': True, 'new_external_spending': False,
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'continuation_source_sha256': hashlib.sha256((ROOT/'research/interacting_scaling_20260915/screen_continue.py').read_bytes()).hexdigest()})


def envelope(case):
    protocol = json.loads(PROTOCOL.read_text())['main_proof_if_screen_succeeds']
    if (case/'prepared').exists(): raise ValueError('Do not alter an existing coefficient map')
    path = case/'design.json'; design = json.loads(path.read_text())
    design.update(max_Gram_entries=protocol['Gram_entry_cap'], max_coefficient_rows=protocol['coefficient_row_cap'])
    path.write_text(json.dumps(design, indent=2)+'\n')
    dump(case/'posthoc_capacity.json', {'protocol_sha256': hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        'not_the_frozen_capacity_fallback': True, 'supports_and_ideal_equations_unchanged': True})


def run():
    protocol = json.loads(PROTOCOL.read_text())
    if hashlib.sha256((OUT/'cold/h12_heldout/result.json').read_bytes()).hexdigest() != protocol['frozen_result_sha256']:
        raise ValueError('The original frozen result changed')
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != protocol['driver_sha256']:
        raise ValueError('Declared diagnostic driver changed')
    if hashlib.sha256((ROOT/'research/interacting_scaling_20260915/screen_continue.py').read_bytes()).hexdigest() != protocol['continuation_source_sha256']:
        raise ValueError('Declared continuation changed')
    frozen = json.loads((OUT/'heldout_protocol.json').read_text())
    for name, digest in frozen['source_sha256'].items():
        if hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/name).read_bytes()).hexdigest() != digest:
            raise ValueError('The preserved mathematical procedure changed')
    from research.interacting_scaling_20260915.screen_continue import run as continuation
    start = time.monotonic(); deadline = start+protocol['maximum_additional_seconds']
    source = OUT/'cases'/protocol['source_screen']
    attempts, selected = [], None
    for number, seconds in enumerate(protocol['same_family_continuation_slices_seconds']):
        name = 'h12_posthoc_screen_'+str(number)
        record = continuation(source, name, seconds, deadline)
        attempts.append(record)
        if not record['execution']['completed']: break
        source = OUT/'cases'/name
        if record['exact_screen']['meets_full_interval_lower_threshold']:
            selected = source; break
    result = {'magnetic_continuations': attempts, 'selected_screen': str(selected) if selected else None,
        'frozen_H12_result_unchanged': True, 'fresh_complete_cold_calculation': False,
        'original_failed_frozen_run_seconds_additional': json.loads((OUT/'cold/h12_heldout/result.json').read_text())['wall_seconds']}
    if selected:
        case = OUT/'cases/h12_posthoc_collective4'
        spec = protocol['main_proof_if_screen_succeeds']
        stages = [('initialize', 20, [STD, '-B', '-S', '-m', PREFIX+'cases', case.name,
            str(OUT/'cases/h12_heldout_base'), '--widths']+list(map(str, spec['widths']))+['--collective-widths']+list(map(str, spec['collective_widths']))),
            ('attach', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'attach', str(case), str(selected)]),
            ('capacity', 20, [STD, '-B', '-S', '-m', PREFIX+'posthoc_h12', 'envelope', str(case)]),
            ('prepare', 900, [NUM, '-B', '-m', PREFIX+'prepare', str(case)]),
            ('solve', spec['solve_seconds']+45, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', str(spec['solve_seconds'])]),
            ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'complete', str(OUT/'models/h12_heldout'), str(case), 'solve', str(case/'exact')])]
        result['main_execution'] = execute(case.name, stages, deadline)
        if result['main_execution']['completed']:
            result['interval'] = json.loads((case/'original_interval.json').read_text())
    result['additional_measured_continuation_seconds'] = time.monotonic()-start
    dump(OUT/'posthoc_h12_result.json', result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('declare', 'run', 'envelope')); p.add_argument('case', nargs='?', type=Path)
    a = p.parse_args()
    declare() if a.action == 'declare' else run() if a.action == 'run' else envelope(a.case.resolve())
