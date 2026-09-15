"""Exact rational orbital localization proposal and finite-sector diagnostics."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_orbital_rotation import rotate, car_check
from experiments.marginal_symbolic import decode, encode
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import select_reference, spin_states, spin_blocks


def matrix_product(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b))] for i in range(len(a))]


def rational_orthogonal(target, denominator=10000):
    """Approximate numerical Givens angles, preserving orthogonality exactly."""
    import numpy as np
    if type(denominator) is not int or not 1 <= denominator <= 100000:
        raise ValueError('Bounded positive half-angle denominator required')
    a = np.asarray(target, dtype=float)
    if (a.ndim != 2 or a.shape[0] != a.shape[1] or not 1 <= len(a) <= 16
            or not np.isfinite(a).all() or np.max(np.abs(a.T@a-np.eye(len(a)))) > 1e-7):
        raise ValueError('Finite numerical orthogonal target of bounded dimension required')
    n = len(a); a = a.copy(); steps = []
    for column in range(n-1):
        for row in range(n-1, column, -1):
            p, q = row-1, row
            x, y = a[p, column], a[q, column]
            norm = float(np.hypot(x, y))
            if not norm: continue
            c, s = float(x/norm), float(y/norm)
            first, second = a[p].copy(), a[q].copy()
            a[p], a[q] = c*first+s*second, -s*first+c*second
            # Either chart keeps the half-angle parameter bounded near pi.
            h = F(s/(1+c) if c >= 0 else s/(1-c)).limit_denominator(denominator)
            cosine = (1-h*h)/(1+h*h) if c >= 0 else (h*h-1)/(1+h*h)
            sine = 2*h/(1+h*h)
            steps.append((p, q, cosine, sine))
    diagonal = [F(1 if a[i,i] >= 0 else -1) for i in range(n)]
    u = [[F(int(i==j)) for j in range(n)] for i in range(n)]
    for p, q, c, s in steps:
        inverse = [[F(int(i==j)) for j in range(n)] for i in range(n)]
        inverse[p][p] = inverse[q][q] = c
        inverse[p][q], inverse[q][p] = -s, s
        u = matrix_product(u, inverse)
    u = [[u[i][j]*diagonal[j] for j in range(n)] for i in range(n)]
    if not car_check(u): raise ValueError('Rational construction lost orthogonality')
    return u, {'givens_count': len(steps), 'half_angle_denominator_cap': denominator,
               'maximum_target_entry_error': float(np.max(np.abs(np.array(u,dtype=float)-np.asarray(target)))),
               'maximum_entry_numerator_bits': max(x.numerator.bit_length() for row in u for x in row),
               'maximum_entry_denominator_bits': max(x.denominator.bit_length() for row in u for x in row)}


def replay_rotation(certificate):
    if certificate.get('kind') != 'rational_spatial_rotation_v1':
        raise ValueError('Unsupported orbital rotation certificate')
    base = {k: certificate[k] for k in ('modes','particles','hamiltonian')}
    original = SpinZeroOracle(base)
    u = [[F(x) for x in row] for row in certificate['rotation']]
    if len(u)*2 != original.modes or not car_check(u):
        raise ValueError('Rotation dimension or CAR identity failed')
    transformed = rotate(original.h, u)
    declared = decode(certificate['rotated_hamiltonian'], original.modes, 4)
    if transformed != declared: raise ValueError('Transformed Hamiltonian coefficients do not match')
    inverse = [list(row) for row in zip(*u)]
    if rotate(transformed, inverse) != original.h: raise ValueError('Exact inverse recovery failed')
    rotated = SpinZeroOracle(dict(base,hamiltonian=encode(transformed)))
    return {'car_valid': True, 'inverse_exact': True, 'spin_symmetry_valid': True,
            'modes': original.modes, 'particles': original.particles,
            'original_terms': len(original.h), 'rotated_terms': len(rotated.h),
            'scope': 'Exact number-preserving real spatial orbital rotation applied equally to both spins; CAR, forward coefficient equality, inverse recovery and spin symmetry checked. This is a unitary change of the Hamiltonian representation, not an energy or complement-gap certificate.'}


def sector_diagnostic(data, seed_policy='minimum_diagonal'):
    import numpy as np
    oracle = SpinZeroOracle(data); states = spin_states(oracle)
    if seed_policy not in ('minimum_diagonal','default'):
        raise ValueError('Unknown retained-space seed policy')
    seed = min(states, key=lambda s: (oracle.action(s).get(s,F(0)),s)) if seed_policy=='minimum_diagonal' else None
    p, _, history = select_reference(oracle, 32, seed=seed)
    blocks = spin_blocks(oracle, p)
    full = np.array([[float(oracle.action(t).get(s,F(0))) for t in states] for s in states])
    spectrum = np.linalg.eigvalsh(full)
    physical = []; ceilings = {str(k): [] for k in (2,3,4)}
    for group in blocks:
        a = np.array([[float(oracle.action(t).get(s,F(0))) for t in group] for s in group])
        physical.append(float(np.linalg.eigvalsh(a)[0]))
        for k in (2,3,4):
            comparison = -np.abs(a)/(k-1);np.fill_diagonal(comparison,np.diag(a))
            ceilings[str(k)].append(float(np.linalg.eigvalsh(comparison)[0]))
    q_floor = min(physical)
    return spectrum, {'seed_policy': seed_policy, 'seed': p[0], 'p_states': p,
        'retained_upper': history[-1]['upper'], 'block_dimensions': [len(b) for b in blocks],
        'ground_energy_numerical': float(spectrum[0]), 'physical_q_floor_numerical': q_floor,
        'physical_q_gap_numerical': q_floor-float(spectrum[0]),
        'comparison_family_ceilings_numerical': {k: min(v) for k,v in ceilings.items()},
        'scope': 'Explicit spin-sector numerical diagnostic with fresh retained-space selection. Comparison-family values are not certified physical lower bounds; no scalability or exact cone-feasibility claim.'}


def construct(source, target_source, output, denominator=10000):
    source, target_source, output = Path(source), Path(target_source), Path(output)
    if output.exists(): raise ValueError('Preserve previous localized basis export')
    data = json.loads(source.read_text()); base = {k:data[k] for k in ('modes','particles','hamiltonian')}
    proposal = json.loads(target_source.read_text())
    if proposal.get('source_sha256') != hashlib.sha256(source.read_bytes()).hexdigest():
        raise ValueError('Localization target must bind the source Hamiltonian file')
    u, fitting = rational_orthogonal(proposal['target_rotation'], denominator)
    original = SpinZeroOracle(base); transformed = rotate(original.h,u)
    certificate = dict(base,kind='rational_spatial_rotation_v1',
                       rotation=[[rational_text(x) for x in row] for row in u],
                       rotated_hamiltonian=encode(transformed))
    receipt = replay_rotation(certificate)
    rotated = dict(base,hamiltonian=encode(transformed))
    e0, before = sector_diagnostic(base); e1, after = sector_diagnostic(rotated)
    _, default_after = sector_diagnostic(rotated, 'default')
    import numpy as np
    diagnostic = {'source':str(source),'target_source':str(target_source),'fitting':fitting,
                  'spectrum_max_abs_error':float(np.max(np.abs(e0-e1))),
                  'original':before,'localized':after,'localized_default_seed_control':default_after,
                  'scope':'Exact rational approximation to the numerically aligned localized orbital target, with exact rotation replay and finite numerical gap diagnostics. No localized energy certificate or general scaling result.'}
    output.mkdir(parents=True)
    for name, value in [('rotation_certificate.json',certificate),('rotation_receipt.json',receipt),
                        ('hamiltonian.json',rotated),('diagnostic.json',diagnostic)]:
        (output/name).write_text(json.dumps(value,indent=2)+'\n')
    return diagnostic


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--source');parser.add_argument('--target-source')
    parser.add_argument('--output');parser.add_argument('--verify');parser.add_argument('--denominator',type=int,default=10000)
    args=parser.parse_args()
    if args.verify: print(json.dumps(replay_rotation(json.loads(Path(args.verify).read_text())),indent=2))
    elif args.source and args.target_source and args.output:
        print(json.dumps(construct(args.source,args.target_source,args.output,args.denominator),indent=2))
    else: parser.error('Specify --verify or --source, --target-source and --output')
