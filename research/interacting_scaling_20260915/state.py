"""Fresh globally charge-correct product start for local-orbital MPS discovery."""
import argparse
import json
from pathlib import Path
from research.interacting_scaling_20260915.budget import dump
from research.molecular_collective_20260913.core import digest
from research.nvidia_followup_20260915.strict_replay import require_supported_sector


def seed(data):
    require_supported_sector(data)
    spatial, count = data['modes']//2, data['particles']//2
    orders = [list(range(parity, spatial, 2))+list(range(1-parity, spatial, 2)) for parity in (0, 1)]
    occupied = [set(order[:count]) for order in orders]
    charges, tensors, current = [[[0, 0]]], [], [0, 0]
    for i in range(data['modes']):
        bit = int(i//2 in occupied[i%2])
        tensors.append([[0, bit, 0, 1]])
        current = current.copy()
        current[i%2] += bit
        charges.append([current])
    return {'kind': 'integer_charge_mps_v1', 'fixture_sha256': digest(data),
        'modes': data['modes'], 'particles': data['particles'], 'spin_counts': [count, count],
        'denominator': 1, 'bond_charges': charges, 'tensors': tensors}


def run(case, bond=32, sweeps=6, initial=None):
    data = json.loads((case/'fixture.json').read_text())
    if initial is None:
        initial = case/'product_seed.json'
        dump(initial, seed(data))
    from research.correlated_pair_20260913.mps_spatial import run as discover
    folder = case/'mps'
    folder.mkdir(exist_ok=False)
    discover(data, [data['particles']//2]*2, folder, bond=bond, sweeps=sweeps, initial=str(initial))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('--bond', type=int, default=32)
    p.add_argument('--sweeps', type=int, default=6)
    p.add_argument('--initial', type=Path)
    a = p.parse_args()
    run(a.case.resolve(), a.bond, a.sweeps, a.initial)
