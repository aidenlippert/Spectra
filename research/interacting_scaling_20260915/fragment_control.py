"""Exact charge-complete local control for the declared zero-coupling model.

Only each two-spatial-orbital fragment is enumerated. All fragment charges and
magnetic sectors are included; global number is imposed by min-plus convolution.
This control is not a constructor input for the interacting proof.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from experiments.marginal_symbolic import canonical, decode, encode
from research.compact_response_20260913.closure import check_factor
from research.transfer_followup_20260915.enumerated_baseline import blocks, quotient, shifted
from research.interacting_scaling_20260915.budget import dump


def fragments(case):
    data = json.loads((case/'fixture.json').read_text())
    receipt = json.loads((case/'coupling.json').read_text())
    if F(receipt['coupling']) != 0: raise ValueError('This control only applies at exactly zero coupling')
    partition = receipt['partition']
    if sorted(i for part in partition for i in part) != list(range(data['modes']//2)) or any(len(part) > 2 for part in partition):
        raise ValueError('Disjoint complete fragments of at most two spatial orbitals required')
    owner = {i: k for k, part in enumerate(partition) for i in part}
    h = decode(data['hamiltonian'], data['modes'], 4)
    terms = [{} for _ in partition]
    for word, value in h.items():
        touched = {owner[i//2] for c, i in word}
        if len(touched) > 1: raise ValueError('Interfragment coupling cannot enter the disconnected control')
        k = next(iter(touched)) if touched else 0
        index = {p: q for q, p in enumerate(partition[k])}
        local = tuple((c, 2*index[i//2]+i%2) for c, i in word)
        terms[k][local] = value
    reconstructed = {}
    for part, poly in zip(partition, terms):
        for word, value in poly.items():
            global_word = tuple((c, 2*part[i//2]+i%2) for c, i in word)
            reconstructed[global_word] = reconstructed.get(global_word, F(0))+value
    if canonical(reconstructed) != h: raise ValueError('Exact local Hamiltonian sum failed')
    return data, partition, terms


def combine(energies, total):
    table = {0: (F(0), [])}
    for local in energies:
        next_table = {}
        for count, (value, choices) in table.items():
            for charge, bound in enumerate(local):
                if count+charge > total: continue
                candidate = value+bound, choices+[charge]
                if count+charge not in next_table or candidate[0] < next_table[count+charge][0]:
                    next_table[count+charge] = candidate
        table = next_table
    if total not in table: raise ValueError('Global particle number cannot be attained')
    return table[total]


def construct(case, output):
    import numpy as np
    data, partition, polynomials = fragments(case)
    records = []
    for part, h in zip(partition, polynomials):
        rows = []
        for charge in range(2*len(part)+1):
            local = {'modes': 2*len(part), 'particles': charge, 'hamiltonian': encode(h)}
            groups, denominator, _ = blocks(local)
            for alpha, states, H in groups:
                A = np.array(H, dtype=float)/denominator
                values, vectors = np.linalg.eigh(A)
                lower = F(round(float(values[0])*10**12), 10**12)-F(1, 10**8)
                fd = 10**12
                C = np.linalg.cholesky(A-(float(lower)+1e-10)*np.eye(len(H)))
                factor = [[round(float(C[i, j])*fd) for j in range(i+1)] for i in range(len(H))]
                vector = [round(float(v)*fd) for v in vectors[:, 0]]
                rows.append({'charge': charge, 'alpha': alpha, 'lower_Ha': str(lower),
                    'factor_denominator': fd, 'factor': factor, 'upper_vector': vector})
        records.append(rows)
    dump(output/'certificate.json', {'partition': partition, 'blocks': records,
        'all_local_charges_and_magnetic_sectors_included': True,
        'used_by_interacting_state_or_lower_discovery': False})


def check(case, output):
    import sys
    if any(name in sys.modules for name in ('numpy', 'scipy', 'pyscf')):
        raise ValueError('Numerical imports are not allowed in standalone fragment acceptance')
    data, partition, polynomials = fragments(case)
    certificate = json.loads((output/'certificate.json').read_text())
    if certificate['partition'] != partition or len(certificate['blocks']) != len(partition):
        raise ValueError('Missing or changed fragment')
    lower_curves, upper_curves, labels = [], [], 0
    for part, h, supplied in zip(partition, polynomials, certificate['blocks']):
        expected, lower, upper = [], [], []
        for charge in range(2*len(part)+1):
            local = {'modes': 2*len(part), 'particles': charge, 'hamiltonian': encode(h)}
            groups, denominator, _ = blocks(local)
            lows, ups = [], []
            for alpha, states, H in groups:
                expected.append((charge, alpha)); labels += len(states)
                row = next((r for r in supplied if (r['charge'], r['alpha']) == (charge, alpha)), None)
                if row is None: raise ValueError('An allowed fragment charge/spin block is missing')
                L = F(row['lower_Ha'])
                K, kd = shifted(H, denominator, L)
                check_factor(K, kd, F(0), row['factor'], row['factor_denominator'])
                U = quotient(H, denominator, row['upper_vector'])
                if L > U: raise ValueError('Invalid local endpoints')
                lows.append(L); ups.append(U)
            lower.append(min(lows)); upper.append(min(ups))
        if [(r['charge'], r['alpha']) for r in supplied] != expected:
            raise ValueError('Charge/spin coverage is not exact and unique')
        lower_curves.append(lower); upper_curves.append(upper)
    L, lower_charges = combine(lower_curves, data['particles'])
    U, upper_charges = combine(upper_curves, data['particles'])
    # A common exact supporting slope supplies an additional composable bound.
    # Piecewise-linear maxima occur at a breakpoint of one fragment's envelope.
    slopes = {F(0)}
    for curve in lower_curves:
        slopes.update((curve[q]-curve[p])/(q-p) for p in range(len(curve)) for q in range(p+1, len(curve)))
    supporting, mu = max((sum(min(e-slope*q for q, e in enumerate(curve)) for curve in lower_curves)+slope*data['particles'], slope) for slope in slopes)
    if supporting > L: raise ValueError('Supporting bound exceeded the exact charge-convolution lower')
    dump(output/'replay.json', {'lower_Ha': str(L), 'upper_Ha': str(U), 'width_mHa': float(1000*(U-L)),
        'target_met': U-L <= F(1, 625), 'valid_on': 'Entire global fixed-N sector of the exact disconnected model',
        'local_labels_enumerated_including_all_charges': labels, 'global_determinant_labels_enumerated': 0,
        'global_charge_assignments_enumerated': 0, 'charge_combination': 'Exact min-plus convolution over total particle count',
        'lower_minimizing_charges': lower_charges, 'upper_product_state_charges': upper_charges,
        'supporting_chemical_potential_Ha': str(mu), 'supporting_lower_Ha': str(supporting),
        'supporting_interval_width_mHa': float(1000*(U-supporting)),
        'local_lower_curves_Ha': [[str(v) for v in c] for c in lower_curves],
        'additive_local_approximation_losses_are_retained': True,
        'nonzero_interfragment_coupling_is_refused': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('construct', 'check'))
    p.add_argument('case', type=Path); p.add_argument('output', type=Path); a = p.parse_args()
    (construct if a.action == 'construct' else check)(a.case.resolve(), a.output.resolve())
