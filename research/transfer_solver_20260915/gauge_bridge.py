"""Exact orbital-sign unitary and coefficient-norm bridge between fixtures."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import time


def polynomial(data):
    return {tuple(map(tuple, t['word'])): F(t['coefficient']) for t in data['hamiltonian']}


def bridge(old, new):
    if (old['modes'], old['particles']) != (new['modes'], new['particles']):
        raise ValueError('Gauge comparison requires matching particle/orbital counts')
    a, b = polynomial(old), polynomial(new)
    pivots = {}
    compatible = True
    for w in sorted(a.keys() & b.keys()):
        mask = 0
        for creation, orbital in w:
            mask ^= 1 << (orbital//2)
        rhs = int(a[w]*b[w] < 0)
        while mask:
            pivot = (mask & -mask).bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = mask, rhs
                break
            row, value = pivots[pivot]
            mask ^= row
            rhs ^= value
        if not mask and rhs:
            compatible = False
            break
    signs = 0
    if compatible:
        for pivot in sorted(pivots, reverse=True):
            mask, rhs = pivots[pivot]
            if rhs ^ ((mask & signs).bit_count() % 2):
                signs |= 1 << pivot
    transformed = {w: c*((-1)**sum((signs >> (orbital//2)) & 1 for creation, orbital in w)) for w, c in b.items()}
    residual = sum(abs(a.get(w, F(0))-transformed.get(w, F(0))) for w in a.keys() | b.keys())
    raw = sum(abs(a.get(w, F(0))-b.get(w, F(0))) for w in a.keys() | b.keys())
    return {'paired_spatial_orbital_sign_bits': signs, 'sign_constraints_consistent': compatible,
            'unitary': 'Product of (-1) to the total alpha+beta occupation of each marked spatial orbital',
            'raw_coefficient_l1_difference_Ha': str(raw), 'aligned_coefficient_l1_bound_Ha': str(residual),
            'exact_unitary_equivalence': residual == 0,
            'consequence': 'Any certified [L,U] for the new fixture gives [L-delta,U+delta] for the old fixture on the same N sector.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('old', type=Path)
    p.add_argument('new', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--interval', type=Path)
    a = p.parse_args()
    started = time.monotonic()
    record = bridge(json.loads(a.old.read_text()), json.loads(a.new.read_text()))
    record.update(old_sha256=hashlib.sha256(a.old.read_bytes()).hexdigest(),
                  new_sha256=hashlib.sha256(a.new.read_bytes()).hexdigest())
    if a.interval:
        interval = json.loads(a.interval.read_text())
        delta = F(record['aligned_coefficient_l1_bound_Ha'])
        L, U = F(interval['lower_Ha'])-delta, F(interval['upper_Ha'])+delta
        record.update(original_fixture_lower_Ha=str(L), original_fixture_upper_Ha=str(U),
                      original_fixture_width_mHa=float(1000*(U-L)), target_met=U-L <= F(1, 625))
    record['seconds'] = time.monotonic()-started
    with a.output.open('x') as stream:
        json.dump(record, stream, indent=2)
    print(json.dumps(record, indent=2))
