"""Exact, bounded audit of the collective grammar on frozen molecular inputs.

No energy is solved here. Residual norm bounds are sufficient error budgets,
not lower bounds on the best possible approximation error. The tensor-rank
witness concerns full-Fock-space orbital support, before fixed-N identities.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import decode, mono, add, scale, product, hermitian
from experiments.marginal_general_schur import norm_bound
from experiments.marginal_transfer_verify import apply_word

ROOT = Path(__file__).resolve().parents[2]


def density_word(i, j):
    i, j = sorted((i, j))
    return ((1, i), (1, j), (0, i), (0, j))


def median(values):
    values = sorted(values)
    return values[len(values) // 2]


def residual(h, modes, active):
    active = set(active)
    bath = set(range(modes)) - active
    bb = list(combinations(sorted(bath), 2))
    ab = [(i, j) for i in sorted(active) for j in sorted(bath)]
    # n_i n_j = -a_i† a_j† a_i a_j in the repository's CAR convention.
    jb = median([-h.get(density_word(i, j), F(0)) for i, j in bb])
    jab = median([-h.get(density_word(i, j), F(0)) for i, j in ab])
    kept = {w: c for w, c in h.items()
            if len(w) < 4 or all(i in active for _, i in w)}
    for coefficient, pairs in ((jb, bb), (jab, ab)):
        for i, j in pairs:
            kept = add(kept, mono(density_word(i, j), -coefficient))
    delta = add(h, scale(kept, -1))
    if not hermitian(delta):
        raise AssertionError('Residual lost Hermiticity')
    return delta, jb, jab


def one_leg_rank(h, modes):
    """Exact rank and a nonzero minor of the antisymmetric quartic tensor."""
    rows = [{} for _ in range(modes)]
    for w, c in h.items():
        if len(w) != 4:
            continue
        if tuple(flag for flag, _ in w) != (1, 1, 0, 0):
            raise ValueError('Number-conserving quartic tensor required')
        p, q, r, s = (i for _, i in w)
        for i, j, k, l, sign in ((p, q, r, s, 1), (q, p, r, s, -1),
                                (p, q, s, r, -1), (q, p, s, r, 1)):
            rows[i][(j, k, l)] = c * sign
    echelon = []
    selected_rows = []
    for i, row in enumerate(rows):
        current = dict(row)
        for pivot, previous in echelon:
            factor = current.get(pivot, F(0))
            if factor:
                for column, value in previous.items():
                    value = current.get(column, F(0)) - factor * value
                    if value:
                        current[column] = value
                    else:
                        current.pop(column, None)
        if current:
            pivot = min(current)
            factor = current[pivot]
            echelon.append((pivot, {c: v / factor for c, v in current.items()}))
            selected_rows.append(i)
    columns = [pivot for pivot, _ in echelon]
    minor = [[rows[i].get(c, F(0)) for c in columns] for i in selected_rows]
    determinant = det(minor)
    if not determinant and minor:
        raise AssertionError('Reported tensor-rank minor is singular')
    return {'rank': len(echelon), 'rows': selected_rows,
            'columns': [list(c) for c in columns], 'minor_determinant': str(determinant),
            'minor': [[str(c) for c in row] for row in minor],
            'scope': 'Full-Fock quartic tensor support under orbital rotations only; no fixed-N quotient, collective elimination, approximation, or no-go for general compression.'}


def det(matrix):
    a = [[F(v) for v in row] for row in matrix]
    value = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        value *= a[j][j]
        for i in range(j + 1, len(a)):
            factor = a[i][j] / a[j][j]
            for k in range(j + 1, len(a)):
                a[i][k] -= factor * a[j][k]
    return value


def residual_witness(delta, modes, particles):
    """Find a nonzero exact off-diagonal fixed-N matrix element (small oracle).

    Discovery explicitly visits the small fixture sector. Rechecking the
    resulting pair requires only applying each residual word to one state.
    """
    from research.collective_interference_20260913.collective import states
    basis=states(modes,particles);entries={};actions=0
    for w,c in delta.items():
        if tuple(i for f,i in w if f)==tuple(i for f,i in w if not f):continue
        for s in basis:
            actions+=1;target=apply_word(w,s)
            if target:
                t,sign=target;key=(t,s);entries[key]=entries.get(key,F(0))+c*sign
    if not entries:raise ValueError('No off-diagonal residual witness')
    (bra,ket),value=max(entries.items(),key=lambda item:(abs(item[1]),item[0]))
    witness={'bra':bra,'ket':ket,'matrix_element':str(value)}
    check_residual_witness(delta,modes,particles,witness)
    return {**witness,'operator_norm_lower_bound_Ha':str(abs(value)),
            'discovery_full_sector_dimension':len(basis),'discovery_word_state_actions':actions,
            'scope':'One exact matrix element of the selected residual, invariant under adding a scalar. This is not an energy error or a no-go for other patterns/collective certificates.'}


def check_residual_witness(delta,modes,particles,witness):
    bra=witness['bra'];ket=witness['ket']
    if any(type(s) is not int or not 0<=s<2**modes or s.bit_count()!=particles for s in (bra,ket)) or bra==ket:
        raise ValueError('Invalid fixed-N residual witness states')
    value=F(0)
    for w,c in delta.items():
        target=apply_word(w,ket)
        if target is not None and target[0]==bra:value+=c*target[1]
    if not value or value!=F(witness['matrix_element']):raise ValueError('Residual matrix-element witness failed')
    return abs(value)


def audit(path):
    start = time.monotonic()
    raw = path.read_bytes()
    data = json.loads(raw)
    m = data['modes']
    if m not in (8, 12):
        raise ValueError('Audit is budgeted for the eight/twelve-mode fixtures')
    h = decode(data['hamiltonian'], m, 4)
    rows = []
    for active in combinations(range(m), 4):
        delta, jb, jab = residual(h, m, active)
        eta = norm_bound(delta, 'hermitian_pairs')
        rows.append((eta, active, delta, jb, jab))
    eta, active, delta, jb, jab = min(rows, key=lambda r: (r[0], r[1]))
    offdiagonal = {w: c for w, c in delta.items()
                   if tuple(i for flag, i in w if flag) != tuple(i for flag, i in w if not flag)}
    return {'fixture': str(path.relative_to(ROOT)), 'fixture_sha256': hashlib.sha256(raw).hexdigest(),
            'modes': m, 'particles': data['particles'], 'active_subsets_examined': len(rows),
            'best_active_modes': list(active), 'collective_bath_density_coefficient': str(jb),
            'collective_active_bath_density_coefficient': str(jab),
            'residual_paired_norm_bound_Ha': str(eta), 'residual_raw_l1_Ha': str(norm_bound(delta)),
            'conservative_interval_widening_2eta_Ha': str(2 * eta),
            'widening_over_1_6mHa_budget': str(2 * eta / F(16, 10000)),
            'residual_quartic_terms': len(delta), 'offdiagonal_residual_terms': len(offdiagonal),
            'offdiagonal_paired_norm_bound_Ha': str(norm_bound(offdiagonal, 'hermitian_pairs')),
            'selected_residual_exact_witness': residual_witness(delta,m,data['particles']),
            'one_leg_rank': one_leg_rank(h, m),
            'wall_seconds': time.monotonic() - start,
            'scope': 'Best sufficient coefficient-norm budget over four existing orbitals and two uniform density coefficients. Not a computed energy interval, not an optimal spectral-norm error, no orbital optimization; arbitrary quadratic terms were retained without claiming they fit the two-band grammar.'}


def run(out):
    from research.collective_interference_20260913.collective import model
    from research.collective_interference_20260913.validation import full_polynomial
    start = time.monotonic()
    molecular = [audit(ROOT / f'results/certificate_scaling/active_space_ladder/{name}/fixture.json')
                 for name in ('h4', 'h6')]
    # Full tensor rank does not preclude our own collective fixed-N elimination.
    control, p = full_polynomial(model(2))
    control_rank = one_leg_rank(control, p['physical_modes'])
    # Rank-one coupling with distinct energies generates all B Krylov directions.
    b = 8
    energies = [F(4) + F(j, 2 * (b - 1)) for j in range(b)]
    krylov = [[e ** k for k in range(b)] for e in energies]
    vandermonde = F(1)
    for i, j in combinations(range(b), 2):
        vandermonde *= energies[j] - energies[i]
    if det(krylov) != vandermonde or not vandermonde:
        raise AssertionError('Krylov countercontrol failed')
    result = {'molecular': molecular, 'collectively_eliminated_control_rank': control_rank,
              'dispersed_rank_one_krylov_control': {'bath_modes': b, 'coupling_rank': 1,
                'krylov_rank': b, 'exact_determinant': str(vandermonde)},
              'wall_seconds': time.monotonic() - start}
    out.write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.out)
    for row in result['molecular']:
        print(row['modes'], row['best_active_modes'],
              float(F(row['residual_paired_norm_bound_Ha'])), row['one_leg_rank']['rank'])
    print('seconds', result['wall_seconds'])
