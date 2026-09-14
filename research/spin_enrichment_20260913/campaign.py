"""Continue the existing exact-certified subspace under a whole-run budget."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

from research.spin_completion_20260913.campaign import ROOT, save
from research.molecular_collective_20260913.core import digest

TARGET = F(16, 10000)
MAX_DIRECTIONS = 256
SOURCE = ROOT/'results/spin_completion_20260913/campaign/adaptive/round_4_batch_2'


def validate_source(model, cert, span):
    if cert.get('kind') != 'spin_completion_subspace_v1' or cert.get('fixture_sha256') != digest(model.data) or cert.get('tail_sha256') != digest(model.tail):
        raise ValueError('Inherited certificate input binding failed')
    if not isinstance(span, list) or not span:
        raise ValueError('Empty inherited subspace')
    from research.spin_completion_20260913.discovery import DIRECTION_DENOMINATOR
    groups, _ = model.groups(span)
    if len(groups) != len(cert['anti_blocks']):
        raise ValueError('Inherited direction/certificate block mismatch')
    for group, block in zip(groups, cert['anti_blocks']):
        if (block['generators'] != group['generators'] or block['directions'] != group['directions']
                or block['direction_denominator'] != DIRECTION_DENOMINATOR):
            raise ValueError('Inherited directions do not match the certified subspace')


def accepted_in_budget(row, budget):
    return row.get('accepted') is True and row['complete_at_seconds'] <= budget


def choose_best(rows, upper, budget, inherited_lower=None):
    good = [r for r in rows if accepted_in_budget(r, budget)
        and (inherited_lower is None or F(r['lower_Ha']) >= inherited_lower)]
    if not good:
        return None
    sufficient = [r for r in good if F(0) <= upper-F(r['lower_Ha']) <= TARGET]
    if sufficient:
        return min(sufficient, key=lambda r: (r['directions'], -F(r['lower_Ha'])))['case']
    return max(good, key=lambda r: F(r['lower_Ha']))['case']


def arm(name, out, source, budget, start):
    from research.certificate_scaling.streaming_reference_upper import upper
    from research.spin_completion_20260913.core import replay
    from research.spin_completion_20260913.discovery import Model
    out.mkdir(parents=True, exist_ok=False)
    prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    cert = json.loads((source/'certificate.json').read_text()); span = json.loads((source/'span.json').read_text())
    model = Model(data, tail); validate_source(model, cert, span)
    inherited = replay(data, tail, cert)
    reference_path = ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json'
    reference = json.loads(reference_path.read_text()); U, reference_receipt = upper(data, reference['independent_upper'])
    if U != F(reference['upper']):
        raise ValueError('Reference upper did not reproduce')
    summary = {'arm': name, 'budget_seconds': budget, 'source': str(source.resolve()),
        'source_hashes': {file: hashlib.sha256((source/file).read_bytes()).hexdigest() for file in ('certificate.json', 'span.json', 'dual_proposal.json')},
        'source_replay': inherited, 'source_directions': len(span), 'upper_Ha': str(U), 'reference_replay': reference_receipt,
        'model_construction': model.construction, 'trials': [], 'rounds': [], 'best_case': None,
        'complete': False, 'pending_trial': None}
    save(out/'summary.json', summary)
    current = {'case': 'inherited', 'lower_Ha': inherited['original_lower_Ha'], 'directions': len(span)}
    dual = json.loads((source/'dual_proposal.json').read_text())
    def remaining():
        return budget-(time.monotonic()-start)
    def trial(label, proposed, separate=False):
        before = time.monotonic()
        summary['pending_trial'] = {'case': label, 'directions': len(proposed), 'started_at_seconds': before-start}
        save(out/'summary.json', summary)
        try:
            receipt, proposal = model.solve(proposed, out/label, remaining(), separate=separate, solver_cap=40.)
            lower = F(receipt['accepted']['original_lower_Ha'])
            if lower > U:
                raise ValueError('Inconsistent certified interval')
            row = {'case': label, 'accepted': True, 'directions': len(proposed), 'lower_Ha': str(lower),
                'width_mHa': float(1000*(U-lower)), 'certificate_bytes': receipt['certificate_bytes'],
                'trial_wall_seconds': time.monotonic()-before, 'complete_at_seconds': time.monotonic()-start}
        except Exception as error:
            row = {'case': label, 'accepted': False, 'directions': len(proposed), 'error': repr(error),
                'trial_wall_seconds': time.monotonic()-before, 'complete_at_seconds': time.monotonic()-start}
            proposal = None
        summary['trials'].append(row); summary['pending_trial'] = None
        summary['best_case'] = choose_best(summary['trials'], U, budget, F(inherited['original_lower_Ha']) if name == 'continue' else None)
        summary['elapsed_seconds'] = time.monotonic()-start
        save(out/'summary.json', summary); print(json.dumps(row), flush=True)
        return row, proposal
    if name == 'separate':
        trial('separate_selected', span, separate=True); summary['stop'] = 'ablation_completed'
    elif name == 'continue':
        for rnd in range(1, 13):
            if U-F(current['lower_Ha']) <= TARGET:
                summary['stop'] = 'target_reached'; break
            if remaining() < 60:
                summary['stop'] = 'time_reserve'; break
            candidates, pricing = model.price(dual, span)
            save(out/f'pricing_{rnd}.json', {'pricing': pricing, 'candidates': candidates})
            round_record = {'round': rnd, 'from_case': current['case'], 'pricing': pricing, 'choices': []}
            choices = []; seen = set(); reached = False
            for count in (1, 2):
                extra = [e for group in candidates for e in group[:count]]; proposed = span+extra
                key = json.dumps(proposed, sort_keys=True)
                if not extra or key in seen or len(proposed) > MAX_DIRECTIONS:
                    continue
                seen.add(key)
                if any(sum(e['group'] == gid for e in proposed) > len(frame['polynomials']) for gid, frame in enumerate(model.frames)):
                    continue
                if remaining() < 60:
                    round_record['unrun_bundle'] = count; break
                row, proposed_dual = trial(f'round_{rnd}_batch_{count}', proposed)
                if accepted_in_budget(row, budget):
                    gain = F(row['lower_Ha'])-F(current['lower_Ha'])
                    score = float(gain)/(row['trial_wall_seconds']+pricing['seconds'])
                    choice = {'case': row['case'], 'gain_mHa': float(1000*gain), 'gain_Ha_per_second': score}
                    round_record['choices'].append(choice)
                    if gain >= F(1, 10**6):
                        choices.append((score, proposed, row, proposed_dual))
                    if U-F(row['lower_Ha']) <= TARGET:
                        reached = True; break
            if choices:
                _, span, current, dual = max(choices, key=lambda item: item[0])
                round_record['selected'] = current['case']; summary['continuation_case'] = current['case']
            else:
                round_record['selected'] = None
            summary['rounds'].append(round_record); summary['elapsed_seconds'] = time.monotonic()-start
            if reached:
                summary['stop'] = 'target_reached'
            elif not choices:
                summary['stop'] = 'no_accepted_gain_or_budget'
            save(out/'summary.json', summary)
            if reached or not choices:
                break
        else:
            summary['stop'] = 'round_budget'
    else:
        raise ValueError('Unknown continuation arm')
    summary.update({'complete': True, 'elapsed_seconds': time.monotonic()-start})
    save(out/'summary.json', summary)


def run(out):
    out.mkdir(parents=True, exist_ok=False); start = time.monotonic(); rows = []; source = SOURCE
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    for name, budget in [('continue', 360.), ('separate', 90.)]:
        before = time.monotonic(); failure = None
        cmd = [sys.executable, '-m', 'research.spin_enrichment_20260913.campaign', '--arm', name,
            '--out', str(out/name), '--source', str(source), '--budget', str(budget), '--started', str(before)]
        with (out/f'{name}.log').open('w') as log:
            try:
                proc = subprocess.run(cmd, cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT, timeout=budget)
                if proc.returncode:
                    failure = f'exit_{proc.returncode}'
            except subprocess.TimeoutExpired:
                failure = 'end_to_end_wall_timeout'
        row = {'arm': name, 'budget_seconds': budget, 'wall_seconds': time.monotonic()-before, 'failure': failure}
        path = out/name/'summary.json'
        if path.exists():
            summary = json.loads(path.read_text())
            row.update({'complete': summary['complete'], 'best_case': summary['best_case'],
                'completed_valid_trials': sum(accepted_in_budget(t, budget) for t in summary['trials'])})
            if name == 'continue' and summary['best_case']:
                source = (out/name/summary['best_case']).resolve()
        rows.append(row); save(out/'summary.json', {'arms': rows, 'wall_seconds': time.monotonic()-start})
        print(json.dumps(row), flush=True)
        if name == 'continue' and source == SOURCE:
            break


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--arm', choices=['continue', 'separate']); parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--budget', type=float, default=360.); parser.add_argument('--started', type=float)
    args = parser.parse_args()
    if args.arm:
        arm(args.arm, args.out, args.source, args.budget, args.started if args.started is not None else time.monotonic())
    else:
        run(args.out)
