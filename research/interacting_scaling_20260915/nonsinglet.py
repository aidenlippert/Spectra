"""Direct magnetic-sector screen, with no global cubic dictionary at level zero."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
from research.interacting_scaling_20260915.budget import OUT, dump


def initialize(source, name, level=0):
    from research.interacting_scaling_20260915.dictionary import windows
    if level not in (0, 1, 2): raise ValueError('Unknown magnetic adaptation level')
    case = OUT/'cases'/name
    case.mkdir(parents=True, exist_ok=False)
    for item in ('fixture.json', 'upper.json'):
        shutil.copyfile(source/item, case/item)
    spatial = json.loads((source/'fixture.json').read_text())['modes']//2
    clusters = [c for width in (2, 3) if width <= spatial for c in windows(spatial, width)] if level else []
    dump(case/'design.json', {'clusters': clusters, 'collective_pairs': level == 2, 'complete': False,
        'magnetization': 1, 'number_body': 2 if level else 1, 'max_Gram_entries': 2000000,
        'max_coefficient_rows': 180000, 'global_particle_number_only': True})
    dump(case/'input_dependencies.json', {'source': str(source),
        'guide': 'exact uniform M_S=1 trace; no MPS or cubic moment input',
        'upper_is_stopping_threshold_only': True})
    return case


def check(case, proposal, output):
    from research.collective_completion_20260914.spin_screen import check_sector
    data = json.loads((case/'fixture.json').read_text())
    certificate = json.loads((case/proposal/'export/certificate.json').read_text())
    if certificate.get('magnetization') != 1 or certificate.get('singlet'):
        raise ValueError('A complete M_S=1 sector proof is required')
    result = check_sector(data, certificate)
    threshold = F(json.loads((case/'upper.json').read_text())['upper_Ha'])-F(1, 625)
    result['meets_full_interval_lower_threshold'] = F(result['lower'])-F(result['original_H_spin_defect_Ha']) >= threshold
    result['not_a_full_fixed_N_bound_by_itself'] = True
    dump(output, result)
    print(json.dumps({'lower_float': result['lower_float'], 'threshold_met': result['meets_full_interval_lower_threshold']}), flush=True)


def attach(source, magnetic):
    """Recheck the candidate before making it a complete-pipeline dependency."""
    from research.collective_completion_20260914.spin_screen import check_sector
    from research.molecular_collective_20260913.core import digest
    data = json.loads((source/'fixture.json').read_text())
    if digest(data) != digest(json.loads((magnetic/'fixture.json').read_text())):
        raise ValueError('Magnetic screen belongs to a different input')
    path = magnetic/'solve/export/certificate.json'
    certificate = json.loads(path.read_text())
    if certificate.get('magnetization') != 1 or certificate.get('singlet'):
        raise ValueError('Only a nonsinglet screen may be attached')
    result = check_sector(data, certificate)
    upper = F(json.loads((source/'upper.json').read_text())['upper_Ha'])
    if upper-(F(result['lower'])-F(result['original_H_spin_defect_Ha'])) > F(1, 625):
        raise ValueError('The direct magnetic screen is valid but too weak for the declared interval target')
    if (source/'nonsinglet.json').exists():
        raise ValueError('Preserve the existing magnetic proof; attach to a new case')
    shutil.copyfile(path, source/'nonsinglet.json')
    dump(source/'nonsinglet_dependency.json', {'source': str(magnetic),
        'certificate': str(path), 'exact_replay': result,
        'construction': json.loads((magnetic/'design.json').read_text()),
        'global_full_cubic_map_or_state_input': False})


def construct(source, deadline):
    """Try nested direct screens, charging every failed level to the caller."""
    from research.interacting_scaling_20260915.pipeline import execute, PREFIX
    from research.interacting_scaling_20260915.budget import STD, NUM
    if (source/'nonsinglet.json').exists():
        raise ValueError('Fresh magnetic construction requires a case without an existing screen')
    outcomes, success = [], False
    for level in range(3):
        case = OUT/'cases'/f'{source.name}_magnetic_level{level}'
        stages = [('initialize', 20, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'initialize', str(source), case.name, '--level', str(level)]),
            ('prepare', 300, [NUM, '-B', '-m', PREFIX+'prepare', str(case)]),
            ('solve', 195, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', '150']),
            ('check', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'check', str(case), 'solve', str(case/'exact.json')])]
        result = execute(case.name, stages, deadline)
        item = {'level': level, 'case': str(case), 'execution': result}
        sufficient = result['completed'] and json.loads((case/'exact.json').read_text())['meets_full_interval_lower_threshold']
        item['threshold_met'] = bool(sufficient)
        outcomes.append(item)
        if sufficient:
            installed = execute(case.name+'_selected', [('attach', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'attach', str(source), str(case)])], deadline)
            success = installed['completed']
            item['attachment'] = installed
            break
        if not result['completed']: break
    receipt = {'completed': success, 'levels': outcomes, 'global_cubic_map_constructed': False}
    dump(source/'magnetic_adaptation.json', receipt)
    return receipt


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='action', required=True)
    i = sub.add_parser('initialize'); i.add_argument('source', type=Path); i.add_argument('name')
    i.add_argument('--level', type=int, default=0)
    c = sub.add_parser('check'); c.add_argument('case', type=Path); c.add_argument('proposal'); c.add_argument('output', type=Path)
    a = sub.add_parser('attach'); a.add_argument('source', type=Path); a.add_argument('magnetic', type=Path)
    a = p.parse_args()
    if a.action == 'initialize': print(initialize(a.source.resolve(), a.name, a.level))
    elif a.action == 'check': check(a.case.resolve(), a.proposal, a.output.resolve())
    else: attach(a.source.resolve(), a.magnetic.resolve())
