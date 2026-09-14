"""Bounded, nonaccepting dual-mixture proposal for the nine-shape family."""
from pathlib import Path
from fractions import Fraction as F
import json
import argparse
import sys
import time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from joint_profile_numeric import build
from joint_family_limit_numeric import solve_exact
from experiments.marginal_charged_projectors import charged_vectors

BASE = ROOT/'results/marginal_graded_hubbard8/joint_projector/signed_density'
OUT = BASE/'symmetric_diagonals_9/extra_31'


def main():
    global OUT
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,default=OUT);parser.add_argument('--range-two',action='store_true');parser.add_argument('--free-range-two-profile',action='store_true');parser.add_argument('--quadratic-charge',action='store_true');parser.add_argument('--charge-square-pairs',action='store_true')
    args=parser.parse_args();OUT=args.directory
    started = time.monotonic()
    c = json.loads((OUT/'profile_joint_r1_2_certificate.json').read_text())
    candidates = json.loads((BASE/'symmetry_diagonal_candidates.json').read_text())['candidates']
    proposal_rows=json.loads((OUT/'profile_proposals.json').read_text())['rows']
    indices=proposal_rows[0].get('shape_indices',list(range(8))+[31])
    if len(indices)!=9 or len(set(indices))!=9:raise ValueError('This bounded dual probe requires nine distinct shapes')
    if any(F(c['target'][k])!=v for k,v in [('U',F(4)),('t',F(1)),('V',F(1,2))]):raise ValueError('This dual probe is for U4,t1,V1/2')
    ys = [{int(s): v for s, v in candidates[i]['diagonal'].items()} for i in indices]
    if ('W' in c['target'])!=args.range_two:raise ValueError('Explicit range-two mode must match the target')
    W=F(c['target']['W']) if args.range_two else None
    if args.free_range_two_profile and not args.range_two:raise ValueError('Free profile requires explicit range-two mode')
    if args.quadratic_charge and not args.free_range_two_profile:raise ValueError('Quadratic family requires free range-two profiles')
    if bool(c.get('quadratic_charge_telescope')) != args.quadratic_charge:raise ValueError('Quadratic mode must match energy certificate')
    quadratic_terms=None
    if args.quadratic_charge:
        from experiments.marginal_quadratic_charge_telescope import coefficients,local_value
        quadratic_terms=coefficients(c['quadratic_charge_telescope'])
        if set(quadratic_terms)!={'0,3'}:raise ValueError('Numeric probe requires the independent quadratic direction')
    if args.charge_square_pairs and not args.quadratic_charge:raise ValueError('Square pairs require full quadratic mode')
    if bool(c.get('charge_square_pair_telescope'))!=args.charge_square_pairs:raise ValueError('Square-pair mode must match energy certificate')
    square_terms=None
    if args.charge_square_pairs:
        from experiments.marginal_charge_square_pairs import LABELS as SQUARE_LABELS, coefficients as square_coefficients, local_value as square_value
        square_terms=square_coefficients(c['charge_square_pair_telescope'])
    range_profile=None
    if args.range_two:
        from experiments.marginal_range_two_density import local_profile,diagonal_value
        range_profile=local_profile(6,W,c['local_window']['range_two_density_profile'])
        if not args.free_range_two_profile and range_profile!=local_profile(6,W):raise ValueError('Fixed probe requires the uniform local profile')
    diagonal_count=9+int(args.free_range_two_profile)+int(args.quadratic_charge)+4*int(args.charge_square_pairs)
    direction_count=6+diagonal_count
    equality_count=1+direction_count
    row_count=equality_count+2
    local = c['local_window']
    x = [F(local[name][i]) for name, i in [('onsite_profile', 0), ('onsite_profile', 1), ('hopping_profile', 0), ('hopping_profile', 1), ('density_profile', 0), ('density_profile', 1)]]
    half = {int(s): v for s, v in c['vector'].items()}
    hn = sum(v*v for v in half.values())
    charged, cn = charged_vectors(c['joint']['vector'])
    ratio = F(c['joint']['ratio'])
    alpha, beta = F(c['penalty']), F(c['joint']['penalty'])
    y = {int(s): F(v) for s, v in c['telescoping_diagonal'].items()}
    original = build(x,W=W,range_profile=range_profile,quadratic_terms=quadratic_terms,square_terms=square_terms)
    perturbed = []
    for j in range(6):
        point = x.copy()
        point[j] += F(1, 100)
        perturbed.append(build(point,W=W,range_profile=range_profile,quadratic_terms=quadratic_terms,square_terms=square_terms))
    rows, states, energies, exact_data = [], [], [], {}
    rng = np.random.default_rng(20260912)
    for key, (_, k, cols) in original.items():
        scale = np.sqrt([sum(v*v for v in col.values()) for col in cols])
        denominator = scale[:, None]*scale[None, :]
        a = np.array(k, float)/denominator
        deriv_exact = [[[100*(v-k[i][j]) for j, v in enumerate(row)] for i, row in enumerate(p[key][1])] for p in perturbed]
        deriv = np.array(deriv_exact, float)/denominator
        ts = np.array([[shape.get(next(iter(col)) & 1023, 0)-shape.get(next(iter(col)) >> 2, 0) for col in cols] for shape in ys])
        if args.free_range_two_profile:ts=np.vstack([ts,[int(diagonal_value(next(iter(col)),[1,-1,-1,1])) for col in cols]])
        if args.quadratic_charge:ts=np.vstack([ts,[int(local_value(next(iter(col)),{'0,3':F(1)})) for col in cols]])
        if args.charge_square_pairs:ts=np.vstack([ts,[[int(square_value(next(iter(col)),{label:F(1)})) for col in cols] for label in SQUARE_LABELS]])
        correction = np.array([float(y.get(next(iter(col)) & 1023, 0)-y.get(next(iter(col)) >> 2, 0)) for col in cols])
        def projector(v, norm):
            w = np.array([sum(b*v.get(s, 0) for s, b in col.items()) for col in cols], float)/scale/np.sqrt(float(norm))
            return np.outer(w, w)
        ph = projector(half, hn)
        q = ph+float(ratio)*sum((projector(v, cn) for v in charged), np.zeros_like(a))
        m = a+float(alpha)*ph+float(beta)*q+np.diag(correction)
        value = eigh(m, subset_by_index=[0, 0], eigvals_only=True)[0]
        if value > float(F(c['penalized_lower']))+1e-4:
            continue
        exact_data[key] = (k, deriv_exact, cols, ts)
        directions = [np.zeros(direction_count)]+[sign*np.eye(direction_count)[j]*.002 for j in range(direction_count) for sign in (-1, 1)]
        directions += [rng.normal(size=direction_count)*.002 for _ in range(16)]
        for delta in directions:
            matrix = m+np.einsum('i,ijk->jk', delta[:6], deriv)+np.diag(delta[6:]@ts)
            _, vectors = eigh(matrix, subset_by_index=[0, 0])
            v = vectors[:, 0]/scale
            coeff = [round(float(z/max(abs(v)))*10**8) for z in v]
            vector = {str(s): z*b for z, col in zip(coeff, cols) if z for s, b in col.items()}
            norm = sum(v*v for v in vector.values())
            w = np.array(coeff, float)*scale/np.sqrt(float(norm))
            d = np.array([w@b@w for b in deriv])
            row = np.r_[1, d, ts@(w*w), w@ph@w, w@q@w]
            if any(np.max(abs(row-old)) < 1e-11 for old in rows):
                continue
            rows.append(row)
            states.append((vector, norm, coeff, key))
            energies.append(float(w@a@w-d@np.array(x, float)))
    th = F(c['projector_sum_ceiling'])/c['windows']
    tj = F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
    arr = np.array(rows).T
    rhs = [F(1)]+[F(0)]*direction_count+[th, tj]
    result = linprog(energies, A_eq=arr[:equality_count], b_eq=np.array(rhs[:equality_count], float), A_ub=arr[equality_count:], b_ub=np.array(rhs[equality_count:], float), bounds=(0, None), method='highs', options={'dual_feasibility_tolerance': 1e-9, 'primal_feasibility_tolerance': 1e-9})
    diagnostic = {'accepted': False, 'states': len(states), 'active_sectors': list(exact_data), 'status': result.message, 'seconds': time.monotonic()-started}
    if result.success:
        values = list(result.x)+list(np.array(rhs[equality_count:], float)-arr[equality_count:]@result.x)
        support = [i for i, v in enumerate(values) if v > 1e-10]
        diagnostic.update(numerical_family_upper=float(result.fun)/5, basis_size=len(support), support=support)
        if len(support) == row_count:
            exact_rows, exact_energies = [], []
            for index in support:
                if index >= len(states):
                    exact_rows.append([F(int(j == equality_count+index-len(states))) for j in range(row_count)])
                    exact_energies.append(F(0))
                    continue
                vector, norm, coeff, key = states[index]
                k, ds, cols, ts = exact_data[key]
                def expectation(matrix):
                    return sum((F(a*b)*matrix[i][j] for i, a in enumerate(coeff) if a for j, b in enumerate(coeff) if b and matrix[i][j]), F(0))/norm
                d = [expectation(matrix) for matrix in ds]
                t = [F(sum(int(ts[j, i])*z*z*sum(b*b for b in cols[i].values()) for i, z in enumerate(coeff)), norm) for j in range(diagonal_count)]
                def fidelity(source, source_norm):
                    return F(sum(v*source.get(int(s), 0) for s, v in vector.items())**2, norm*source_norm)
                ph = fidelity(half, hn)
                q = ph+ratio*sum((fidelity(v, cn) for v in charged), F(0))
                exact_rows.append([F(1)]+d+t+[ph, q])
                exact_energies.append(expectation(k)-sum(a*b for a, b in zip(x, d)))
            weights = solve_exact([[row[i] for row in exact_rows] for i in range(row_count)], rhs)
            if min(weights) < 0:
                raise ValueError('Proposed exact basis has negative weights')
            upper = sum(w*e for w, e in zip(weights, exact_energies))/5
            proposal = {'kind': 'joint_diagonal_family_limit_proposal_v1', 'half_vector': c['vector'], 'charged_vector': c['joint']['vector'], 'ratio': str(ratio), 'theta_half': str(th), 'theta_joint': str(tj), 'diagonal_shapes': [{str(s): str(v) for s, v in shape.items()} for shape in ys], 'mixture': [{'weight': str(w), 'vector': states[i][0]} for i, w in zip(support, weights) if i < len(states)], 'proposed_periodic_family_upper': str(upper), 'scope': 'Untrusted exact-basis proposal. No production verifier yet accepts these extra moment constraints.'}
            if args.range_two:proposal.update(kind='joint_diagonal_range2_family_limit_proposal_v1',W=str(W),range_two_density_profile=c['local_window']['range_two_density_profile'])
            if args.free_range_two_profile:proposal['kind']='joint_diagonal_range2_free_profile_proposal_v1'
            if args.quadratic_charge:proposal['kind']='joint_diagonal_quadratic_charge_family_proposal_v1'
            if args.charge_square_pairs:proposal['kind']='joint_diagonal_charge_square_pair_family_proposal_v1'
            (OUT/'diagonal_family_limit_proposal.json').write_text(json.dumps(proposal, indent=2)+'\n')
            diagnostic.update(rational_family_upper=str(upper), rational_family_upper_float=float(upper), mixture_sources=len(proposal['mixture']), seconds=time.monotonic()-started)
    (OUT/'diagonal_family_limit_diagnostic.json').write_text(json.dumps(diagnostic, indent=2)+'\n')
    print(json.dumps(diagnostic), flush=True)


if __name__ == '__main__':
    main()
