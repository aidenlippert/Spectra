"""One declared residual-phase diagnostic after the fixed H12 continuations.

This changes an optimizer phase, never a proof family or acceptance equation.
Both diagnostic parent clocks together are capped at 2400 additional seconds.
The failed frozen cold run and its discovery costs remain additional.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.pipeline import execute, PREFIX

PROTOCOL = OUT/'residual_h12_protocol.json'


def declare():
    dump(PROTOCOL, {'declared_UTC': datetime.now(timezone.utc).isoformat(),
        'purpose': 'Test the existing mu=0.03 residual phase when the automatic switch has not fired',
        'gate': 'Both fixed continuations completed without an accepted magnetic screen; last phase remains mu=2',
        'source_case': 'h12_posthoc_screen_1', 'refinement_solve_seconds': 300,
        'maximum_combined_diagnostic_parent_seconds': 2400,
        'earlier_frozen_cold_seconds_additional': True,
        'not_a_frozen_transfer_or_fresh_complete_calculation': True,
        'operator_family_and_exact_acceptance_unchanged': True,
        'main_proof_if_screen_passes': 'Same local/collective width-4 recipe as posthoc_h12_protocol.json',
        'new_external_spending': False,
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'prior_protocol_sha256': hashlib.sha256((OUT/'posthoc_h12_protocol.json').read_bytes()).hexdigest()})


def run():
    protocol = json.loads(PROTOCOL.read_text())
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != protocol['driver_sha256']:
        raise ValueError('Declared refinement driver changed')
    if hashlib.sha256((OUT/'posthoc_h12_protocol.json').read_bytes()).hexdigest() != protocol['prior_protocol_sha256']:
        raise ValueError('Prior diagnostic protocol changed')
    frozen = json.loads((OUT/'heldout_protocol.json').read_text())
    for name, digest in frozen['source_sha256'].items():
        if hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen constructor changed')
    prior = json.loads((OUT/'posthoc_h12_result.json').read_text())
    source = OUT/'cases'/protocol['source_case']
    result = {'earlier_diagnostic_parent_seconds': prior['additional_measured_continuation_seconds'],
        'original_failed_cold_seconds_additional': prior['original_failed_frozen_run_seconds_additional'],
        'fresh_complete_cold_calculation': False, 'frozen_H12_result_unchanged': True}
    if prior['selected_screen'] or not (source/'exact.json').exists():
        result['skipped'] = 'A magnetic screen already passed or the required fixed continuations did not complete'
    elif float(json.loads((source/'solve/history.json').read_text())[-1]['mu']) != 2.:
        result['skipped'] = 'The existing residual phase was already exercised'
    else:
        start = time.monotonic()
        available = protocol['maximum_combined_diagnostic_parent_seconds']-result['earlier_diagnostic_parent_seconds']
        deadline = start+max(0, available)
        case = OUT/'cases/h12_posthoc_residual_screen'
        result['screen_execution'] = execute(case.name, [
            ('initialize', 30, [STD, '-B', '-S', '-m', PREFIX+'screen_continue', 'initialize', str(source), case.name]),
            ('solve', 345, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', '300', '--mu', '0.03', '--restart', str(case/'resume.npz')]),
            ('check', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'check', str(case), 'solve', str(case/'exact.json')])], deadline)
        if result['screen_execution']['completed']:
            result['exact_screen'] = json.loads((case/'exact.json').read_text())
            if result['exact_screen']['meets_full_interval_lower_threshold']:
                target = OUT/'cases/h12_posthoc_residual_collective4'
                stages = [
                    ('initialize', 20, [STD, '-B', '-S', '-m', PREFIX+'cases', target.name, str(OUT/'cases/h12_heldout_base'), '--widths', '2', '3', '4', '--collective-widths', '3', '4']),
                    ('attach', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'attach', str(target), str(case)]),
                    ('capacity', 20, [STD, '-B', '-S', '-m', PREFIX+'posthoc_h12', 'envelope', str(target)]),
                    ('prepare', 900, [NUM, '-B', '-m', PREFIX+'prepare', str(target)]),
                    ('solve', 495, [NUM, '-B', '-m', PREFIX+'solve', str(target), 'solve', '--seconds', '450']),
                    ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'complete', str(OUT/'models/h12_heldout'), str(target), 'solve', str(target/'exact')])]
                result['main_execution'] = execute(target.name, stages, deadline)
                if result['main_execution']['completed']:
                    result['interval'] = json.loads((target/'original_interval.json').read_text())
        result['additional_refinement_parent_seconds'] = time.monotonic()-start
        result['combined_diagnostic_parent_seconds'] = result['earlier_diagnostic_parent_seconds']+result['additional_refinement_parent_seconds']
    dump(OUT/'residual_h12_result.json', result)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('declare', 'run'))
    declare() if p.parse_args().action == 'declare' else run()
