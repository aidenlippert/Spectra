"""Keep numerical FCI and exact enumerated comparisons separate from discovery."""
import argparse
from math import comb
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, CHEM, dump
from research.interacting_scaling_20260915.pipeline import execute


def numerical(model):
    return execute(model.name+'_numerical_fci', [('solve', 300, [CHEM, '-B', '-m',
        'research.transfer_solver_20260915.benchmark', str(model)])])


def exact(source, name):
    data = json.loads((source/'fixture.json').read_text())
    if comb(data['modes'], data['particles']) > 2000:
        raise ValueError('Preserve the established small enumerated baseline envelope')
    case = OUT/'models'/name
    case.mkdir(parents=True, exist_ok=False)
    dump(case/'specification.json', json.loads((source/'specification.json').read_text()))
    start = time.monotonic()
    result = execute(name, [('integrals', 300, [CHEM, '-B', '-m', 'research.transfer_solver_20260915.generate', str(case)]),
        ('construct', 180, [NUM, '-B', '-m', 'research.transfer_followup_20260915.enumerated_baseline', 'construct', str(case)]),
        ('check', 180, [STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.enumerated_baseline', 'verify', str(case)])])
    record = {'execution': result, 'fresh_complete_pipeline_seconds': time.monotonic()-start,
        'source_control': str(source), 'matched_original_input': json.loads((case/'fixture.json').read_text()) == data,
        'independent_state_and_lower_discovery': True, 'FCI_or_Spectra_witness_not_an_input': True}
    if result['completed']: record['interval'] = json.loads((case/'enumerated_baseline/replay.json').read_text())
    dump(case/'baseline_result.json', record)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('numerical', 'exact'))
    p.add_argument('model', type=Path); p.add_argument('--name'); a = p.parse_args()
    numerical(a.model.resolve()) if a.action == 'numerical' else exact(a.model.resolve(), a.name)
