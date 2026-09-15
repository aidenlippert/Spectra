"""Gain-per-cost adaptive search and two serial watchdog-bounded controls."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n'); temporary.replace(path)


def arm(name, out, budget, started, adaptive_source=None, diagnostic_source=None):
    from research.spin_completion_20260913.discovery import Model
    start = started if started is not None else time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    model = Model(data, tail)
    summary = {'arm': name, 'budget_seconds': budget, 'model_construction': model.construction,
        'trials': [], 'rounds': [], 'best_case': None, 'complete': False}
    save(out/'summary.json', summary)
    def remaining():
        return budget-(time.monotonic()-start)
    def trial(label, span, separate=False, solver_cap=40.):
        before = time.monotonic()
        try:
            receipt, dual = model.solve(span, out/label, remaining(), separate, solver_cap)
            row = {'case': label, 'accepted': True, 'directions': len(span),
                'lower_Ha': receipt['accepted']['original_lower_Ha'], 'trial_wall_seconds': time.monotonic()-before,
                'complete_at_seconds': time.monotonic()-start, 'certificate_bytes': receipt['certificate_bytes']}
        except Exception as error:
            row = {'case': label, 'accepted': False, 'directions': len(span), 'error': repr(error),
                'trial_wall_seconds': time.monotonic()-before, 'complete_at_seconds': time.monotonic()-start}
            receipt = dual = None
        summary['trials'].append(row)
        good = [r for r in summary['trials'] if r['accepted'] and r['complete_at_seconds'] <= budget]
        summary['best_case'] = max(good, key=lambda r: F(r['lower_Ha']))['case'] if good else None
        summary['elapsed_seconds'] = time.monotonic()-start
        save(out/'summary.json', summary); print(json.dumps(row), flush=True)
        return row, dual
    if name == 'adaptive':
        from research.spin_completion_20260913.diagnostic import check_separator
        from research.molecular_collective_20260913.core import digest
        diagnostic = json.loads((diagnostic_source/'receipt.json').read_text())
        prior_dual = json.loads((ROOT/'results/spin_subspace_20260913/full_dual/witness.json').read_text())
        if diagnostic['fixture_sha256'] != digest(data) or diagnostic['tail_sha256'] != digest(tail) or diagnostic['dual_sha256'] != digest(prior_dual):
            raise ValueError('Diagnostic input mismatch')
        summary['diagnostic_source'] = str(diagnostic_source.resolve())
        span = []; current, dual = trial('baseline', span)
        if not current['accepted']:
            raise RuntimeError('Fresh baseline failed')
        for rnd in range(1, 5):
            if remaining() < 40:
                summary['stop'] = 'time_reserve'; break
            if rnd == 1:
                before = time.monotonic(); candidates = diagnostic['candidates']
                for gid, group in enumerate(candidates):
                    for entry in group:
                        cert = json.loads((diagnostic_source/entry['separator']).read_text())
                        if cert['generators'] != model.frames[gid]['generators'] or cert['vector'] != entry['vector'] or entry['group'] != gid:
                            raise ValueError('Separator/frame ordering mismatch')
                        check_separator(data, tail, prior_dual, cert)
                pricing = {'seconds': time.monotonic()-before, 'source': 'Accepted full-spin counterexample separators; prior discovery charged separately'}
            else:
                candidates, pricing = model.price(dual, span)
            save(out/f'pricing_{rnd}.json', {'pricing': pricing, 'candidates': candidates})
            choices = []; seen = set(); round_record = {'round': rnd, 'pricing': pricing, 'from_case': current['case'], 'choices': []}
            for per_group in (1, 2):
                extra = [e for group in candidates for e in group[:per_group]]
                proposed = span+extra
                key = json.dumps(proposed, sort_keys=True)
                if not extra or key in seen or len(proposed) > 64:
                    continue
                seen.add(key)
                if remaining() < 40:
                    round_record['unrun_bundle'] = per_group; break
                row, proposed_dual = trial(f'round_{rnd}_batch_{per_group}', proposed)
                if row['accepted']:
                    gain = F(row['lower_Ha'])-F(current['lower_Ha'])
                    score = float(gain)/(row['trial_wall_seconds']+pricing['seconds'])
                    choice = {'case': row['case'], 'gain_Ha': str(gain), 'gain_mHa': float(1000*gain), 'gain_Ha_per_second': score}
                    round_record['choices'].append(choice)
                    if gain >= F(1, 10**6):
                        choices.append((score, proposed, row, proposed_dual))
            if not choices:
                round_record['selected'] = None; summary['rounds'].append(round_record)
                summary['stop'] = 'no_accepted_gain_or_budget'; break
            _, span, current, dual = max(choices, key=lambda item: item[0])
            round_record['selected'] = current['case']; summary['rounds'].append(round_record)
            summary['continuation_case'] = current['case']; summary['elapsed_seconds'] = time.monotonic()-start
            save(out/'summary.json', summary)
        else:
            summary['stop'] = 'round_budget'
    elif name == 'full':
        trial('full_frame', model.full_span(), solver_cap=150.)
    elif name == 'separate':
        source = json.loads((adaptive_source/'summary.json').read_text())
        if source['best_case'] is None:
            raise ValueError('No accepted small subspace for ablation')
        source_path = adaptive_source/source['best_case']/'span.json'
        span = json.loads(source_path.read_text())
        absolute_source = source_path.resolve()
        summary['descriptor_source'] = str(absolute_source.relative_to(ROOT) if absolute_source.is_relative_to(ROOT) else absolute_source)
        source_watchdog = json.loads((adaptive_source.parent/'summary.json').read_text())
        summary['descriptor_discovery_cost_seconds'] = next(row['wall_seconds'] for row in source_watchdog['arms'] if row['arm'] == 'adaptive')
        trial('separate_selected', span, separate=True, solver_cap=50.)
    else:
        raise ValueError('Unknown campaign arm')
    summary.update({'complete': True, 'elapsed_seconds': time.monotonic()-start})
    save(out/'summary.json', summary)


def run(out, diagnostic_source):
    out.mkdir(parents=True, exist_ok=False); start = time.monotonic(); rows = []
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    for name, budget in [('adaptive', 240.), ('separate', 90.), ('full', 240.)]:
        before = time.monotonic()
        cmd = [sys.executable, '-m', 'research.spin_completion_20260913.campaign', '--arm', name,
            '--out', str(out/name), '--budget', str(budget), '--started', str(before), '--diagnostic-source', str(diagnostic_source.resolve())]
        if name == 'separate':
            cmd += ['--adaptive-source', str(out/'adaptive')]
        failure = None
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
                'completed_valid_trials': sum(t['accepted'] and t['complete_at_seconds'] <= budget for t in summary['trials'])})
        rows.append(row); save(out/'summary.json', {'arms': rows, 'wall_seconds': time.monotonic()-start})
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--arm', choices=['adaptive', 'separate', 'full']); parser.add_argument('--budget', type=float, default=240.)
    parser.add_argument('--started', type=float); parser.add_argument('--adaptive-source', type=Path)
    parser.add_argument('--diagnostic-source', type=Path, required=True)
    args = parser.parse_args()
    if args.arm:
        arm(args.arm, args.out, args.budget, args.started, args.adaptive_source, args.diagnostic_source)
    else:
        run(args.out, args.diagnostic_source)
