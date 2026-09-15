"""Summarize durable run receipts without treating numerical proposals as bounds."""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
from research.interacting_scaling_20260915.budget import ROOT, OUT, dump


def read(path):
    return json.loads(path.read_text())


def belongs_to_cold(attempt, name):
    """Match declared cold stage names, excluding similarly named later controls."""
    prefix = name+'_'
    if not attempt.startswith(prefix): return False
    suffix = attempt[len(prefix):]
    return bool(re.match(r'(?:input_|state_screen_|base_magnetic_|level\d+_|reference_init_|reference_magnetic_)', suffix)
        or suffix in ('reference_prepare', 'reference_solve', 'reference_replay'))


def family(case):
    prepared = case/'prepared'
    if not (prepared/'frame.json').exists(): return None
    frame = read(prepared/'frame.json')
    result = {'Gram_entries': sum(b['dimension']**2 for b in frame['blocks']),
        'largest_block': max(b['dimension'] for b in frame['blocks']),
        'coefficient_rows': len(frame['rows']), 'coefficient_nonzeros': frame['coefficient_nonzeros'],
        'prepared_bytes': sum(p.stat().st_size for p in prepared.iterdir() if p.is_file()),
        'map_seconds': frame['map_build_seconds'], 'moment_seconds': frame['moment_seconds']}
    if (case/'exact/certificate.json').exists():
        certificate = read(case/'exact/certificate.json')
        result['exported_factor_rows'] = sum(len(block['factor']) for block in certificate['core']['blocks'])
        result['exported_factor_scalar_entries'] = sum(sum(len(row) for row in block['factor']) for block in certificate['core']['blocks'])
        paths = [case/'exact/certificate.json', case/'exact/nonsinglet.json', case/'mps/state.json', case/'upper.json']
        if (case/'rotation.json').exists(): paths.append(case/'rotation.json')
        elif (case/'local_guide.json').exists():
            local = Path(read(case/'local_guide.json')['source'])
            paths.append(local/'rotation.json')
        result['retained_witness_bytes_excluding_shared_input'] = sum(p.stat().st_size for p in paths)
        result['counted_witness_files'] = {str(p.relative_to(ROOT)): p.stat().st_size for p in paths}
        result['witness_size_excludes'] = 'Hamiltonian descriptions, shared checker/source code, and derived replay logs; complete campaign file sizes are in the manifest.'
        state = read(case/'mps/state.json')
        result['MPS_nonzero_integer_entries'] = sum(map(len, state['tensors']))
        result['MPS_maximum_spin_bond'] = max(map(len, state['bond_charges']))
    return result


def run():
    attempts = [read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    unfinished = [d['name'] for d in attempts if d['status'] == 'starting']
    if unfinished: raise ValueError('Complete all measured attempts before final accounting: '+str(unfinished))
    cold = {}
    for path in sorted((OUT/'cold').glob('*/result.json')):
        name = path.parent.name
        result = read(path)
        owned = [d for d in attempts if belongs_to_cold(d['name'], name)]
        child_seconds = sum(d['wall_seconds'] for d in owned)
        if child_seconds > result['wall_seconds']+1:
            raise ValueError('Child stage ledger exceeds its measured serial cold clock: '+name)
        best = result.get('best')
        item = {'target_met': result['target_met'], 'complete_wall_seconds': result['wall_seconds'],
            'sum_child_stage_wall_seconds': child_seconds,
            'peak_child_RSS_bytes': max((d['peak_child_RSS_bytes'] for d in owned), default=0),
            'peak_proof_preparation_RSS_bytes': max((d['peak_child_RSS_bytes'] for d in owned if d['name'].endswith('_prepare')), default=0),
            'attempt_count': len(owned), 'failed_attempts': [d['name'] for d in owned if d['status'] != 'passed'],
            'all_attempts': [d['name'] for d in owned], 'source_result': str(path.relative_to(ROOT))}
        if best:
            interval = best['interval']
            exact_width = F(interval['upper_Ha'])-F(interval['lower_Ha'])
            if exact_width < 0: raise ValueError('An accepted interval has reversed endpoints')
            if bool(exact_width <= F(1, 625)) != bool(interval['target_met']):
                raise ValueError('Exact interval target flag disagrees')
            case = Path(best['case'])
            item.update(width_mHa=float(1000*(F(interval['upper_Ha'])-F(interval['lower_Ha']))),
                accepted_case=str(case.relative_to(ROOT)), family=family(case))
        item['levels'] = [{'level': row['level'], 'case': row['case'], 'family': family(Path(row['case'])),
            'width_mHa': row.get('interval', {}).get('width_mHa')} for row in result['outcomes'] if 'level' in row]
        prepared_cases = sorted({row['case'] for row in result['outcomes']
            if 'case' in row and (Path(row['case'])/'prepared/frame.json').exists()})
        item['sum_main_prepared_Gram_entries_including_rejected_levels'] = sum(
            sum(b['dimension']**2 for b in read(Path(case)/'prepared/frame.json')['blocks'])
            for case in prepared_cases)
        item['prepared_Gram_sum_scope'] = ('Sum over independently prepared main-proof cases; '
            'not peak allocation, stored witness size, runtime, or the separate magnetic-screen work.')
        cold[name] = item
    comparisons = {}
    for size in ('h8', 'h10'):
        direct, reference = size+'_matched_direct', size+'_matched_reference'
        if size+'_matched_reference_extended' in cold: reference = size+'_matched_reference_extended'
        if direct not in cold or reference not in cold: continue
        a, b = cold[direct], cold[reference]
        original_a, original_b = (OUT/'models'/name/'fixture.json' for name in (direct, reference))
        same = read(original_a) == read(original_b)
        row = {'same_exact_fresh_original_input': same, 'both_target_met': a['target_met'] and b['target_met'],
            'one_sample_each_not_a_variance_estimate': True, 'same_host_and_numerical_backends': True,
            'direct_run': direct, 'reference_run': reference}
        if same and row['both_target_met']:
            row['reference_time_divided_by_direct_time'] = b['complete_wall_seconds']/a['complete_wall_seconds']
            row['reference_preparation_peak_divided_by_direct_peak'] = b['peak_proof_preparation_RSS_bytes']/a['peak_proof_preparation_RSS_bytes']
            row['reference_Gram_entries_divided_by_direct'] = b['family']['Gram_entries']/a['family']['Gram_entries']
        comparisons[size] = row
    intervals = {}
    for case in sorted((OUT/'cases').iterdir()):
        path = case/'original_interval.json'
        original = path.exists()
        if not original: path = case/'exact/interval.json'
        if path.exists() and (case/'exact/complete_replay.json').exists():
            d = read(path)
            width = F(d['upper_Ha'])-F(d['lower_Ha'])
            if width < 0: raise ValueError('An accepted case has reversed endpoints')
            intervals[case.name] = {'width_mHa': float(1000*width), 'target_met': width <= F(1, 625),
                'scope': 'Original model with rotation allowance' if original else 'Exact case Hamiltonian',
                'receipt': str(path.relative_to(ROOT))}
    enumerated_agreement = {}
    for path in sorted((OUT/'models').glob('*_enumerated_control/baseline_result.json')):
        control = read(path)
        source = Path(control['source_control']).name
        if source not in cold or not control.get('interval'): continue
        record = read(OUT/'cold'/source/'result.json')
        if not record.get('best'): continue
        matched = control['matched_original_input']
        row = {'same_exact_original_input': matched, 'comparison_receipt': str(path.relative_to(ROOT))}
        if matched:
            a, b = record['best']['interval'], control['interval']
            lower = max(F(a['lower_Ha']), F(b['lower_Ha']))
            upper = min(F(a['upper_Ha']), F(b['upper_Ha']))
            if lower > upper: raise ValueError('Exact independently constructed intervals disagree: '+source)
            row.update(intervals_overlap_exactly=True, overlap_width_mHa=float(1000*(upper-lower)))
        enumerated_agreement[source] = row
    result = {'assembled_UTC': datetime.now(timezone.utc).isoformat(), 'cold_runs': cold,
        'matched_comparisons': comparisons, 'campaign_attempt_count': len(attempts),
        'independent_enumerated_interval_agreement': enumerated_agreement,
        'campaign_child_stage_seconds': sum(d['wall_seconds'] for d in attempts),
        'campaign_failed_attempts': [d for d in attempts if d['status'] != 'passed'],
        'exact_complete_case_intervals': intervals,
        'accuracy_target_misses_are_separate_from_process_failures': [name for name, d in intervals.items() if not d['target_met']],
        'peak_child_RSS_bytes': max(d['peak_child_RSS_bytes'] for d in attempts),
        'ledger_scope': 'Scientific child stages, tests and dependency-install attempts. Parent orchestration clocks are separate; do not add them again. Editorial work and routine file inspection are not timed compute stages.',
        'current_accounting_invocation_excluded_from_its_own_ledger': True,
        'new_cloud_instances': 0, 'new_external_spending_USD': 0,
        'receipt_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((OUT/'runs').glob('*.json'))}}
    dump(OUT/'accounting.json', result)
    print(json.dumps({k: v for k, v in result.items() if k in ('cold_runs', 'matched_comparisons', 'campaign_attempt_count', 'campaign_child_stage_seconds')}, indent=2))


if __name__ == '__main__': run()
