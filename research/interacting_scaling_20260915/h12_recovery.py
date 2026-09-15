"""Correct the post hoc H12 certificate handoff without editing prior receipts."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.pipeline import execute, PREFIX


def run():
    prior = json.loads((OUT/'residual_h12_result.json').read_text())
    if not prior['exact_screen']['meets_full_interval_lower_threshold']:
        raise ValueError('A sufficient exactly checked magnetic screen is required')
    if prior['main_execution']['completed'] or prior['main_execution']['last_stage'] != 'initialize':
        raise ValueError('This recovery concerns only the recorded initialization failure')
    frozen = json.loads((OUT/'heldout_protocol.json').read_text())
    for name, digest in frozen['source_sha256'].items():
        if hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen constructor changed')
    certificate = OUT/'cases/h12_posthoc_residual_screen/solve/export/certificate.json'
    remaining = 2400-prior['combined_diagnostic_parent_seconds']
    case = OUT/'cases/h12_posthoc_main_fixed'
    stages = [
        ('initialize', 20, [STD, '-B', '-S', '-m', PREFIX+'cases', case.name,
            str(OUT/'cases/h12_heldout_base'), '--widths', '2', '3', '4',
            '--collective-widths', '3', '4', '--nonsinglet-source', str(certificate)]),
        ('capacity', 20, [STD, '-B', '-S', '-m', PREFIX+'posthoc_h12', 'envelope', str(case)]),
        ('prepare', 900, [NUM, '-B', '-m', PREFIX+'prepare', str(case)]),
        ('solve', 495, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', '450']),
        ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'complete', str(OUT/'models/h12_heldout'), str(case), 'solve', str(case/'exact')])]
    dump(OUT/'h12_recovery_protocol.json', {'declared_UTC': datetime.now(timezone.utc).isoformat(),
        'reason': 'Pass the newly verified magnetic certificate through the existing --nonsinglet-source argument; the frozen base has none',
        'failed_case_and_original_drivers_preserved': True,
        'new_family_or_acceptance_change': False, 'fresh_cold_calculation': False,
        'earlier_diagnostic_parent_seconds': prior['combined_diagnostic_parent_seconds'],
        'maximum_remaining_parent_seconds': remaining,
        'certificate_source': str(certificate),
        'certificate_sha256': hashlib.sha256(certificate.read_bytes()).hexdigest(),
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'all_failed_frozen_and_diagnostic_costs_additional': True,
        'new_external_spending': False})
    start = time.monotonic()
    result = {'execution': execute(case.name, stages, start+max(0, remaining)),
        'source_case': str(case), 'fresh_complete_cold_calculation': False}
    result['additional_recovery_parent_seconds'] = time.monotonic()-start
    result['combined_diagnostic_parent_seconds'] = prior['combined_diagnostic_parent_seconds']+result['additional_recovery_parent_seconds']
    result['frozen_run_seconds_additional'] = prior['original_failed_cold_seconds_additional']
    if result['execution']['completed']:
        result['interval'] = json.loads((case/'original_interval.json').read_text())
    dump(OUT/'h12_recovery_result.json', result)


if __name__ == '__main__': run()
