"""Run comparison work after the frozen discovery campaign, one process at a time."""
import datetime
import hashlib
import json
import os
import subprocess
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, CHEM, PREFIX


def bounded(name, seconds, command):
    result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', name,
                             '--seconds', str(seconds), '--']+command, cwd=ROOT)
    return {'name': name, 'exit_code': result.returncode}


def run():
    if not (OUT/'validation_set_execution.json').exists():
        raise ValueError('Finish the frozen discovery set before comparison or reference work')
    os.environ['NUMBA_CACHE_DIR'] = str(OUT/'followup_numba_cache')
    sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for folder in ('transfer_solver_20260915', 'transfer_followup_20260915', 'ch2_model_study_20260915')
               for p in (ROOT/'research'/folder).glob('*.py')}
    dump(OUT/'followup_protocol.json', {'frozen_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'source_hashes': sources, 'references_are_post_discovery_only': True,
         'controls': ['numerical FCI for all generated cases', 'independent exact H4 oracle',
                      'Full-sector rigorous enumerated baseline on H4, H6 and water',
                      'H4 and H8 equivalent-cone magnetic controls', 'H8 restricted64 nested-family control', 'exact fragment composition',
                      'CH2 numerical geometry and basis study', 'exact projector component profiling']})
    outcomes = []
    outcomes.append(bounded('h8_cold_fast_projector', 900, [NUM, '-B', '-m', PREFIX+'fast_projector',
                             str(OUT/'cases/h8_cold'), str(OUT/'cases/h8_cold/fast_projector.json')]))
    size = json.loads((OUT/'cases/h10_size/execution.json').read_text())
    if size['status'] != 'completed' and size['last_step'] in ('state', 'upper'):
        result = subprocess.run([STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.adaptive_size', 'run'], cwd=ROOT)
        outcomes.append({'name': 'h10_smaller_state_retry', 'exit_code': result.returncode})
    reference_cases = [OUT/'cases'/name for name in ('h4_control', 'h8_cold', 'h6_asymmetric', 'water_asymmetric', 'h10_size')]
    if (OUT/'adaptive/h10_bond64').exists():
        reference_cases.append(OUT/'adaptive/h10_bond64')
    for case in reference_cases:
        name = case.name
        if (case/'integrals.npz').exists():
            outcomes.append(bounded(name+'_fci', 180, [CHEM, '-B', '-m', PREFIX+'benchmark', str(case)]))
    outcomes.append(bounded('independent_h4_oracle', 180, [STD, '-B', '-S', '-m', PREFIX+'independent_oracle',
        str(OUT/'cases/h4_control'), str(OUT/'independent_h4_oracle.json')]))
    outcomes.append(bounded('fragment_composition', 180, [STD, '-B', '-S', '-m',
        'research.transfer_followup_20260915.fragments']))
    for name in ('h4_control', 'h6_asymmetric', 'water_asymmetric'):
        case = OUT/'cases'/name
        if (case/'fixture.json').exists():
            proposal = bounded(name+'_enumerated_construct', 180, [NUM, '-B', '-m',
                               'research.transfer_followup_20260915.enumerated_baseline', 'construct', str(case)])
            outcomes.append(proposal)
            if proposal['exit_code'] == 0:
                outcomes.append(bounded(name+'_enumerated_verify', 180, [STD, '-B', '-S', '-m',
                                'research.transfer_followup_20260915.enumerated_baseline', 'verify', str(case)]))
    for name in ('h4_control', 'h8_cold'):
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'representation_control',
                                 'run', str(OUT/'cases'/name)], cwd=ROOT)
        outcomes.append({'name': name+'_magnetic_control', 'exit_code': result.returncode})
    result = subprocess.run([STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.restricted_control', 'run'], cwd=ROOT)
    outcomes.append({'name': 'h8_restricted64_control', 'exit_code': result.returncode})
    result = subprocess.run([STD, '-B', '-S', '-m', 'research.ch2_model_study_20260915.study', 'run'], cwd=ROOT)
    outcomes.append({'name': 'ch2_model_study', 'exit_code': result.returncode})
    for name in ('h10_size',):
        case = OUT/'cases'/name
        if (case/'prepared/twirl.npz').exists():
            outcomes.append(bounded(name+'_fast_projector', 900, [NUM, '-B', '-m', PREFIX+'fast_projector',
                                      str(case), str(case/'fast_projector.json')]))
    dump(OUT/'followup_execution.json', outcomes)


if __name__ == '__main__':
    run()
