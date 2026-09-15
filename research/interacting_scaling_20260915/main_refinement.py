"""One final declared same-family residual-phase test for the H12 main proof."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.pipeline import execute, PREFIX

PROTOCOL = OUT/'h12_main_refinement_protocol.json'
SOURCE = OUT/'cases/h12_posthoc_main_fixed'
TARGET = OUT/'cases/h12_main_refined'


def declare():
    dump(PROTOCOL, {'declared_UTC': datetime.now(timezone.utc).isoformat(),
        'source': str(SOURCE), 'reason': 'The unchanged residual phase repaired the magnetic screen; test it once on the unconverged main proof',
        'gate': 'Prior main replay completed, missed 1.6 mHa and never entered mu=0.03',
        'maximum_additional_parent_seconds': 960, 'solve_seconds': 300, 'replay_cap_seconds': 600,
        'additional_to_original_2400_second_diagnostic_envelope': True,
        'not_a_frozen_result_or_matched_resource_comparison': True,
        'old_operator_maps_and_MPS_reused_explicitly': True,
        'all_earlier_failed_attempts_and_source_costs_additional': True,
        'no_new_family_or_acceptance_equation': True, 'new_external_spending': False,
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})


def initialize():
    frame = json.loads((SOURCE/'prepared/frame.json').read_text())
    design = json.loads((SOURCE/'design.json').read_text())
    if frame.get('magnetization') != 0 or design.get('magnetization', 0) != 0:
        raise ValueError('Only an identical singlet-family continuation is supported')
    for field, path in [('fixture_sha256', SOURCE/'fixture.json'), ('state_sha256', SOURCE/'mps/state.json')]:
        if frame[field] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError('Prepared maps or guide moments have a different source')
    if not (SOURCE/'original_interval.json').exists():
        raise ValueError('The prior main proof must complete its independent replay')
    TARGET.mkdir(parents=True, exist_ok=False)
    files = ['fixture.json', 'upper.json', 'nonsinglet.json', 'rotation.json', 'design.json']
    for name in files: shutil.copyfile(SOURCE/name, TARGET/name)
    shutil.copytree(SOURCE/'mps', TARGET/'mps')
    shutil.copytree(SOURCE/'prepared', TARGET/'prepared')
    shutil.copyfile(SOURCE/'solve/checkpoint.npz', TARGET/'resume.npz')
    bound = [SOURCE/name for name in files]+[SOURCE/'mps/state.json', SOURCE/'solve/checkpoint.npz']+sorted((SOURCE/'prepared').iterdir())
    dump(TARGET/'continuation_dependency.json', {'source': str(SOURCE),
        'same_model_sector_state_operators_constraints_and_maps': True,
        'new_coefficient_map_construction': False,
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in bound if p.is_file()},
        'source_discovery_and_all_failed_attempts_additional': True})


def run():
    protocol = json.loads(PROTOCOL.read_text())
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != protocol['driver_sha256']:
        raise ValueError('Declared driver changed')
    frozen = json.loads((OUT/'heldout_protocol.json').read_text())
    for name, digest in frozen['source_sha256'].items():
        if hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen mathematics or acceptance code changed')
    prior = json.loads((OUT/'h12_recovery_result.json').read_text())
    result = {'prior_diagnostic_parent_seconds': prior['combined_diagnostic_parent_seconds'],
        'frozen_cold_seconds_additional': prior['frozen_run_seconds_additional'],
        'fresh_complete_cold_calculation': False}
    if not prior.get('interval') or prior['interval']['target_met']:
        result['skipped'] = 'The prior complete replay failed or already met the target'
    elif float(json.loads((SOURCE/'solve/history.json').read_text())[-1]['mu']) != 2.:
        result['skipped'] = 'The existing residual phase was already exercised'
    else:
        start = time.monotonic()
        stages = [('initialize', 30, [STD, '-B', '-S', '-m', PREFIX+'main_refinement', 'initialize']),
            ('solve', 345, [NUM, '-B', '-m', PREFIX+'solve', str(TARGET), 'solve', '--seconds', '300', '--mu', '0.03', '--restart', str(TARGET/'resume.npz')]),
            ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'complete', str(OUT/'models/h12_heldout'), str(TARGET), 'solve', str(TARGET/'exact')])]
        result['execution'] = execute(TARGET.name, stages, start+protocol['maximum_additional_parent_seconds'])
        result['additional_parent_seconds'] = time.monotonic()-start
        result['all_H12_diagnostic_parent_seconds'] = result['prior_diagnostic_parent_seconds']+result['additional_parent_seconds']
        if result['execution']['completed']:
            result['interval'] = json.loads((TARGET/'original_interval.json').read_text())
    dump(OUT/'h12_main_refinement_result.json', result)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('declare', 'initialize', 'run'))
    {'declare': declare, 'initialize': initialize, 'run': run}[p.parse_args().action]()
