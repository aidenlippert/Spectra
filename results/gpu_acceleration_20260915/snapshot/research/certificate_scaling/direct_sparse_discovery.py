"""Direct bounded sparse squares, optimized jointly with a number multiplier.

Only H,M,N enter discovery. All CAR coefficients, including generated degree-six
terms, enter the residual LP. A frozen H-overlap ranking selects two-word atoms;
this scans a quadratic word-pair pool, and is not a safe omitted-family oracle.
The numerical proposal is never a proof: production Fraction replay is final.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from experiments.marginal_coefficient import dictionaries, dagger, hopping_model
from experiments.marginal_symbolic import (
    add, adj, canonical, decode, encode, hermitian, mono, multiplier_basis,
    number_shift, product, scale, verify, word_product,
)


def square(words, signs):
    return add(*(scale(dict(word_product(dagger(u), v)), F(a*b))
                 for u, a in zip(words, signs) for v, b in zip(words, signs)))


def candidates(h, modes, budget):
    """Generate supports without consulting a solved certificate or wavefunction."""
    atoms, pairs = [], []
    scanned = 0
    for group in dictionaries(modes, 'quadratic'):
        words = group['words']
        for i, u in enumerate(words):
            atoms.append(([u], [1], square([u], [1])))
            for v in words[i+1:]:
                scanned += 1
                cross = add(dict(word_product(dagger(u), v)),
                            dict(word_product(dagger(v), u)))
                score = abs(sum(float(c) * float(h.get(w, 0)) for w, c in cross.items()))
                if score:
                    pairs.append((score, u, v))
    # Deterministic frozen ordering; both signs treated symmetrically.
    atoms.sort(key=lambda a: (-sum(abs(float(h.get(w, 0)*c)) for w, c in a[2].items()), a[0]))
    retained = atoms[:budget]
    # Reserve up to half the slots for off-diagonal structure.
    if pairs and budget:
        retained = atoms[:min(len(atoms), max(1, budget//2))]
        for _, u, v in sorted(pairs, key=lambda a: (-a[0], a[1], a[2])):
            for sign in (1, -1):
                if len(retained) >= budget:
                    break
                retained.append(([u, v], [1, sign], square([u, v], [1, sign])))
            if len(retained) >= budget:
                break
    return retained, {'word_pair_candidates_scanned': scanned,
                      'nonzero_H_overlap_pairs': len(pairs),
                      'monomial_candidates': len(atoms), 'retained_atoms': len(retained)}


def sparse_columns(polys, lookup):
    rr, cc, vv = [], [], []
    for j, poly in enumerate(polys):
        for w, value in poly.items():
            if value:
                rr.append(lookup[w]); cc.append(j); vv.append(float(value))
    return sparse.csc_matrix((vv, (rr, cc)), shape=(len(lookup), len(polys)))


def discover(h, modes, particles, budget=512, ideal_body=1,
             denominator=10**6, time_limit=90):
    start = time.monotonic()
    h = canonical(h)
    if not hermitian(h) or any(len(w)>4 or sum(2*c-1 for c, _ in w) for w in h):
        raise ValueError('Expected Hermitian, number-conserving two-body H')
    if not 0 <= particles <= modes or budget < 0 or ideal_body not in (0, 1, 2):
        raise ValueError('Invalid sector or search budget')
    atoms, counts = candidates(h, modes, budget)
    basis = multiplier_basis(modes, max_body=ideal_body)
    shift = number_shift(modes, particles)
    ideal = [product(shift, p) for p in basis]
    polys = [mono(())] + ideal + [a[2] for a in atoms]
    words = sorted(set(h) | {w for p in polys for w in p}, key=lambda w: (len(w), w))
    lookup = {w: i for i, w in enumerate(words)}
    A = sparse_columns(polys, lookup)
    rhs = np.array([float(h.get(w, 0)) for w in words])
    nfree = 1 + len(ideal)
    nrows = len(words)
    eq = sparse.hstack([A, sparse.eye(nrows), -sparse.eye(nrows)], format='csc')
    # H = b I + (Nhat-N) X + sum_j weight_j B_j^dagger B_j + R.
    # Exact full word basis: c = ||R||_1 - b, without half-row multiplicities.
    objective = np.r_[[-1.], np.zeros(len(ideal)+len(atoms)), np.ones(2*nrows)]
    bounds = [(None, None)]*nfree + [(0, None)]*(len(atoms)+2*nrows)
    built = time.monotonic()
    sol = linprog(objective, A_eq=eq, b_eq=rhs, bounds=bounds,
                  method='highs', options={'time_limit': time_limit})
    solved = time.monotonic()
    if not sol.success:
        raise RuntimeError(f'LP did not converge: {sol.status}: {sol.message}')
    blocks = []
    for j, (ws, signs, _) in enumerate(atoms):
        weight = sol.x[nfree+j]
        integer = int(round(np.sqrt(max(weight, 0))*denominator))
        if integer:
            blocks.append({'name': f'atom-{j}', 'words': ws,
                           'factor': [[integer*s for s in signs]]})
    x = add(*(scale(p, F(int(round(v*10**9)), 10**9))
              for p, v in zip(basis, sol.x[1:nfree])))
    cert = {'modes': modes, 'particles': particles, 'hamiltonian': encode(h),
            'b': str(F(int(round(sol.x[0]*10**9)), 10**9)),
            'number_multiplier': encode(x), 'denominator': denominator, 'blocks': blocks}
    receipt = verify(cert)
    done = time.monotonic()
    residual = rhs - A @ sol.x[:A.shape[1]]
    receipt.update(counts)
    receipt.update({'method': 'direct_H_ranked_two_word_LP', 'atom_budget': budget,
                    'ideal_body': ideal_body, 'ideal_columns': len(ideal),
                    'coefficient_rows': nrows, 'maximum_coefficient_degree': max(map(len, words)),
                    'LP_variables': len(objective), 'LP_nonzeros': eq.nnz,
                    'numeric_proposed_b': float(sol.x[0]),
                    'numeric_physical_residual_l1': float(abs(residual).sum()),
                    'numeric_lower': float(-sol.fun), 'LP_iterations': int(sol.nit),
                    'build_seconds': built-start, 'solve_seconds': solved-built,
                    'export_replay_seconds': done-solved, 'wall_seconds': done-start,
                    'source_factors_used': False, 'source_upper_used': False,
                    'full_Gram_constructed': False, 'full_Fock_space_constructed': False,
                    'omitted_family_optimality_proved': False,
                    'certificate_bytes': len(json.dumps(cert, separators=(',', ':')).encode()),
                    'process_maxrss_native': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    'maxrss_unit': 'bytes' if sys.platform == 'darwin' else 'KiB'})
    return cert, receipt


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fixture', type=Path)
    p.add_argument('--modes', type=int, default=4)
    p.add_argument('--max-atoms', type=int, default=512)
    p.add_argument('--ideal-body', type=int, default=1)
    p.add_argument('--denominator', type=int, default=10**6)
    p.add_argument('--time-limit', type=float, default=90)
    p.add_argument('--outputdir', type=Path, required=True)
    args = p.parse_args()
    if args.fixture:
        source = json.loads(args.fixture.read_text())
        modes, particles = source['modes'], source['particles']
        h = decode(source['hamiltonian'], modes, 4)
    else:
        modes, particles = args.modes, args.modes//2
        h = hopping_model(modes, F(1,5))
    cert, receipt = discover(h, modes, particles, args.max_atoms, args.ideal_body,
                            args.denominator, args.time_limit)
    receipt['H_sha256'] = hashlib.sha256(json.dumps(encode(h), sort_keys=True).encode()).hexdigest()
    args.outputdir.mkdir(parents=True, exist_ok=True)
    (args.outputdir/'certificate.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n')
    (args.outputdir/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
