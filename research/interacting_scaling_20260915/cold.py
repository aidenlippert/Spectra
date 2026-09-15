"""Fresh-integral local-orbital campaign with a frozen nested adaptation rule."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
from research.interacting_scaling_20260915.pipeline import execute, steps_for, PREFIX
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, CHEM, dump


def initialize(name, specification):
    model = OUT/'models'/name
    model.mkdir(parents=True, exist_ok=False)
    dump(model/'specification.json', specification)
    return model


def run(name, specification, levels=((2, 3), (2, 3, 4), (2, 3, 4, 5)), solve_seconds=300, total_seconds=3600, reference_mode=False):
    directory = OUT/'cold'/name
    directory.mkdir(parents=True, exist_ok=False)
    original = OUT/'models'/name
    local = OUT/'cases'/(name+'_base')
    rotation = OUT/'rotations'/name
    source_files = ('cold.py', 'state.py', 'prepare.py', 'dictionary.py', 'spin_patterns.py', 'solve.py',
        'nonsinglet.py', 'warmstart.py', 'pipeline.py', 'cases.py', 'rotation.py', 'budget.py', 'reference.py', 'complete.py')
    source_hashes = {f: hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/f).read_bytes()).hexdigest() for f in source_files}
    dump(directory/'protocol.json', {'specification': specification, 'levels': levels, 'solve_seconds': solve_seconds,
        'complete_pipeline_budget_seconds': total_seconds, 'per_stage_hard_caps_recorded_in_executions': True,
        'new_integrals_and_product_state_discovery': True, 'sources': source_hashes,
        'retained_full_canonical_reference': reference_mode,
        'acceptance': 'Exact complete original-H width <=1.6 mHa; every retained stage charged',
        'continuation': 'Exact checked embedding of earlier local-family proposals; no full-family teacher'})
    start = time.monotonic()
    deadline = start+total_seconds
    model = initialize(name, specification)
    prefix = [('integrals', 300, [CHEM, '-B', '-m', 'research.transfer_solver_20260915.generate', str(model)]),
        ('rotation', 300, [CHEM, '-B', '-m', PREFIX+'rotation', 'construct', str(model), str(rotation)]),
        ('initialize', 20, [STD, '-B', '-S', '-m', PREFIX+'rotation', 'initialize', local.name, str(rotation), '--widths', '2'])]
    outcomes = []
    preparation = execute(name+'_input', prefix, deadline)
    outcomes.append({'input': preparation})
    if not preparation['completed']:
        dump(directory/'result.json', {'outcomes': outcomes, 'target_met': False, 'wall_seconds': time.monotonic()-start})
        return
    # Share fresh state and magnetic-screen construction across nested proof
    # levels within this one cold run; their measured costs remain in its clock.
    common = steps_for(local, True, False, solve_seconds, direct_nonsinglet=True)
    end_label = 'magnetic_initialize'
    common = common[:next(i for i, (label, _, _) in enumerate(common) if label == end_label)]
    preparation = execute(name+'_state_screen', common, deadline)
    outcomes.append({'state_and_screen': preparation})
    best, previous = None, None
    from research.interacting_scaling_20260915.nonsinglet import construct as screen
    if preparation['completed'] and reference_mode:
        case = OUT/'cases'/(name+'_reference')
        result = execute(case.name+'_init', [('reference_initialize', 600,
            [STD, '-B', '-S', '-m', PREFIX+'reference', 'initialize', str(original), str(local), case.name])], deadline)
        outcomes.append({'reference_initialize': result})
        if result['completed']:
            magnetic = screen(case, deadline)
            outcomes.append({'reference_screen': magnetic})
            stages = [('prepare', 900, [NUM, '-B', '-m', PREFIX+'reference', 'prepare', str(case), str(local)]),
                ('solve', min(900, solve_seconds+45), [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', str(solve_seconds)]),
                ('replay', 600, [STD, '-B', '-S', '-m', 'research.nvidia_followup_20260915.strict_replay',
                    str(case), 'solve', str(case/'exact'), '--rotated', str(local)])]
            result = execute(case.name, stages, deadline) if magnetic['completed'] else {'completed': False}
            item = {'case': str(case), 'execution': result}
            if result['completed'] and (case/'exact/complete_replay.json').exists():
                if not json.loads((case/'exact/complete_replay.json').read_text())['all_upper_and_lower_dependencies_rechecked']:
                    raise ValueError('Complete reference replay is required')
                item['interval'] = json.loads((case/'exact/interval.json').read_text())
                best = item
            outcomes.append(item)
    magnetic = screen(local, deadline) if preparation['completed'] and not reference_mode else {'completed': False}
    if not reference_mode: outcomes.append({'direct_screen': magnetic})
    if preparation['completed'] and not reference_mode and magnetic['completed']:
        spatial = json.loads((local/'fixture.json').read_text())['modes']//2
        for level, widths in enumerate(levels):
            if max(widths) > spatial or time.monotonic()-start >= total_seconds: break
            if any(hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/f).read_bytes()).hexdigest() != value for f, value in source_hashes.items()):
                raise ValueError('Frozen constructor changed during this cold run')
            case = OUT/'cases'/f'{name}_level{level}'
            collective = list(map(str, [width for width in widths if width >= 3]))
            command = [STD, '-B', '-S', '-m', PREFIX+'cases', case.name, str(local), '--widths']+list(map(str, widths))
            command += ['--collective-widths']+collective
            stages = [('initialize', 20, command), ('prepare', 900, [NUM, '-B', '-m', PREFIX+'prepare', str(case)])]
            if previous:
                stages.append(('transport', 180, [NUM, '-B', '-m', PREFIX+'warmstart', str(previous), 'solve', str(case)]))
            command = [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', str(solve_seconds)]
            if previous: command += ['--restart', str(case/'nested_start.npz')]
            stages += [('solve', min(900, solve_seconds+45), command),
                ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'complete', str(original), str(case), 'solve', str(case/'exact')])]
            result = execute(case.name, stages, deadline)
            item = {'level': level, 'case': str(case), 'execution': result}
            if result['completed'] and (case/'original_interval.json').exists():
                interval = json.loads((case/'original_interval.json').read_text())
                item['interval'] = interval
                if best is None or interval['width_mHa'] < best['interval']['width_mHa']: best = item
                previous = case
            outcomes.append(item)
            if best and best['interval']['target_met']: break
            if not result['completed']: break
    dump(directory/'result.json', {'outcomes': outcomes, 'best': best,
        'target_met': bool(best and best['interval']['target_met']), 'wall_seconds': time.monotonic()-start,
        'single_fresh_pipeline_on_one_local_host': True, 'all_failed_levels_in_clock': True,
        'integrals_and_state_discovery_included': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('name'); p.add_argument('specification', type=Path)
    p.add_argument('--solve-seconds', type=int, default=300)
    p.add_argument('--reference', action='store_true')
    a = p.parse_args(); run(a.name, json.loads(a.specification.read_text()), solve_seconds=a.solve_seconds, reference_mode=a.reference)
