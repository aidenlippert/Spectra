"""Freeze model choices and the common adaptation rule before new-case results."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from research.interacting_scaling_20260915.budget import ROOT, OUT, dump

PATH = OUT/'heldout_protocol.json'


def freeze():
    geometry = json.loads((ROOT/'results/transfer_solver_20260915/cases/h4_control/specification.json').read_text())['geometry']
    source = ROOT/'research/interacting_scaling_20260915'
    files = ['cold.py', 'reference.py', 'state.py', 'dictionary.py', 'prepare.py', 'spin_patterns.py',
        'solve.py', 'nonsinglet.py', 'pipeline.py', 'warmstart.py', 'rotation.py', 'cases.py', 'budget.py', 'heldout.py', 'complete.py']
    protocol = {'frozen_UTC': datetime.now(timezone.utc).isoformat(),
        'models': {'h12_heldout': {'geometry': [['H', [0., 0., round(1.4*i, 10)]] for i in range(12)], 'basis': 'sto-3g'},
            'water_heldout': {'geometry': [['O', [0., 0., 0.]], ['H', [.77, .03, .57]], ['H', [-.84, .01, .48]]], 'basis': 'sto-3g'},
            'h4_631g_heldout': {'geometry': geometry, 'basis': '6-31g', 'different_model_from_STO3G': True}},
        'levels': [[2, 3], [2, 3, 4], [2, 3, 4, 5]], 'per_level_solve_seconds': 450,
        'complete_pipeline_budget_seconds': 3600, 'local_Gram_entry_cap': 2000000, 'coefficient_row_cap': 180000,
        'magnetic_levels': ['global quadratics', 'quadratics plus local widths 2 and 3',
            'preceding family plus coherent pair-supported cubics'],
        'per_magnetic_level_solve_seconds': 150,
        'source_sha256': {f: hashlib.sha256((source/f).read_bytes()).hexdigest() for f in files},
        'H16_condition': 'Attempt only if H12 meets 1.6 mHa within 1800 s and its measured growth fits the declared local memory/representation envelope.',
        'new_external_spending': False, 'earlier_development_searches_remain_in_campaign_ledger': True}
    dump(PATH, protocol)


def run(name, reference=False):
    protocol = json.loads(PATH.read_text())
    source = ROOT/'research/interacting_scaling_20260915'
    if any(hashlib.sha256((source/f).read_bytes()).hexdigest() != value for f, value in protocol['source_sha256'].items()):
        raise ValueError('The held-out procedure changed after freezing')
    specification = protocol['models'][name]
    from research.interacting_scaling_20260915.cold import run as cold
    cold(name+('_reference' if reference else ''), specification, levels=protocol['levels'],
        solve_seconds=protocol['per_level_solve_seconds'], total_seconds=protocol['complete_pipeline_budget_seconds'], reference_mode=reference)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('freeze', 'run')); p.add_argument('name', nargs='?')
    p.add_argument('--reference', action='store_true'); a = p.parse_args()
    freeze() if a.action == 'freeze' else run(a.name, a.reference)
