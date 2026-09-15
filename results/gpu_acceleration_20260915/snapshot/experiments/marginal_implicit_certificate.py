"""Entire finite-model interval in defect-Dicke coordinates.

Numerics propose coordinates and a lower endpoint. Exact Gram arithmetic
and positivity replay both endpoints without a Fock configuration list.
"""
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import json
import math
import time
from decimal import Decimal, localcontext

from experiments.marginal_defect_dicke import enlarged_workspace, workspace_data, targeted_workspace, combine
from experiments.marginal_enlarged_schur import pivots
from experiments.marginal_general_schur import fixture
from experiments.marginal_schur_transfer import matmul
from experiments.marginal_symbolic import encode, decode


def rational_text(value):
    """Format a checked rational with a local 65,536-bit component budget.

    Small decimal chunks avoid Python's process-wide int/string digit limit;
    the global limit is left unchanged. This is receipt output, not parsing.
    """
    value = F(value)
    def integer_text(number):
        if abs(number).bit_length() > 65536:
            raise ValueError('Rational receipt component exceeds 65536-bit budget')
        sign = '-' if number < 0 else ''
        number = abs(number)
        pieces = []
        while number >= 10**9:
            number, remainder = divmod(number, 10**9)
            pieces.append(f'{remainder:09d}')
        return sign + str(number) + ''.join(reversed(pieces))
    numerator = integer_text(value.numerator)
    return numerator if value.denominator == 1 else numerator + '/' + integer_text(value.denominator)


def orthogonal_coordinates(g):
    n = len(g)
    if not n or any(len(row) != n for row in g) or any(g[i][j] != g[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Expected symmetric square metric')
    lower = [[F(i == j) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        for j in range(i):
            lower[i][j] = (g[i][j] - sum(lower[i][k] * diagonal[k] * lower[j][k] for k in range(j))) / diagonal[j]
        d = g[i][i] - sum(lower[i][k] ** 2 * diagonal[k] for k in range(i))
        if d <= 0:
            raise ValueError('Metric is not positive definite')
        diagonal.append(d)
    transform = [[F(0)] * n for _ in range(n)]
    for j in range(n):
        for i in range(n - 1, -1, -1):
            transform[i][j] = F(i == j) - sum(lower[k][i] * transform[k][j] for k in range(i + 1, n))
    expected = [[diagonal[i] if i == j else F(0) for j in range(n)] for i in range(n)]
    if congruence(g, transform) != expected:
        raise ValueError('Exact orthogonal-coordinate identity failed')
    return transform, diagonal


def congruence(matrix, transform):
    return matmul([list(row) for row in zip(*transform)], matmul(matrix, transform))


def rayleigh(data, coefficients):
    n = data['retained_dimension']
    if type(coefficients) is not list or len(coefficients) != n or any(type(x) is not int for x in coefficients) or not any(coefficients):
        raise ValueError('Expected nonzero integer variational coordinates')
    def quadratic(matrix):
        return sum(coefficients[i] * matrix[i][j] * coefficients[j] for i in range(n) for j in range(n))
    norm = quadratic(data['metric'])
    if norm <= 0:
        raise ValueError('Nonpositive exact variational norm')
    return quadratic(data['projected_h']) / norm


def numerical_coordinates(data, require_leakage=True):
    import numpy as np
    scale = max(data['metric'][i][i] for i in range(len(data['metric'])))
    if scale <= 0:
        raise ValueError('Nonpositive metric scale')
    scaled_metric = [[x / scale for x in row] for row in data['metric']]
    transform, diagonal = orthogonal_coordinates(scaled_metric)
    norms = np.sqrt(np.array(diagonal, dtype=float))
    if not np.all(np.isfinite(norms)) or not np.all(norms > 0):
        raise ValueError('Cannot propose coordinates at this metric scale')
    outer = np.outer(norms, norms)
    a = np.array(congruence([[x / scale for x in row] for row in data['projected_h']], transform), dtype=float) / outer
    k = np.array(congruence([[x / scale for x in row] for row in data['leakage']], transform), dtype=float) / outer if require_leakage else None
    if not np.all(np.isfinite(a)) or (k is not None and not np.all(np.isfinite(k))):
        raise ValueError('Nonfinite projected proposal matrices')
    return transform, norms, a, k


def compress_coordinates(data, coefficients):
    def norm(vector):
        return sum(vector[i] * data['metric'][i][j] * vector[j] for i in range(len(vector)) for j in range(len(vector)))
    original_norm = norm(coefficients)
    if original_norm <= 0:
        raise ValueError('Nonpositive coordinate norm')
    for digits in (12, 18, 24, 36, 48, 72):
        with localcontext() as context:
            context.prec = digits
            rounded = [F(str(Decimal(x.numerator) / Decimal(x.denominator))) for x in coefficients]
        error = norm([a - b for a, b in zip(coefficients, rounded)])
        if error < 0:
            raise ValueError('Negative exact coordinate rounding error')
        if error <= original_norm / 10**24:
            return rounded, {'coordinate_digits': digits,
                             'relative_squared_rounding_error_bound': str(F(1, 10**24)),
                             'relative_squared_rounding_error_float': float(error / original_norm)}
    raise ValueError('Could not compress coordinates within exact physical error budget')


def suggest_upper(data, coordinates=None):
    import numpy as np
    transform, norms, a, _ = coordinates or numerical_coordinates(data, require_leakage=False)
    values, vectors = np.linalg.eigh(a)
    # Significant-digit rounding retains relative precision across metric scales.
    y = [F(format(float(value / norm), '.12g')) for value, norm in zip(vectors[:, 0], norms)]
    x = [sum(t * c for t, c in zip(row, y)) for row in transform]
    x, compression = compress_coordinates(data, x)
    denominator = lcm(*(value.denominator for value in x))
    integers = [int(value * denominator) for value in x]
    divisor = gcd(*integers)
    if not divisor:
        raise ValueError('Rounded variational vector vanished')
    integers = [value // divisor for value in integers]
    upper = rayleigh(data, integers)
    return integers, {'numerical_ritz': float(values[0]), 'upper': str(upper), 'upper_float': float(upper),
                      'maximum_coordinate_bits': max(abs(x).bit_length() for x in integers),
                      'selector': 'Exact Gram LDL with physically checked coordinate compression', **compression}


def suggest_lower(data, coordinates=None):
    import numpy as np
    _, _, a, k = coordinates or numerical_coordinates(data)
    c, eta = float(data['complement_lower']), float(data['norm_bound'])
    def positive(b):
        q = c - eta - b
        return q > 0 and np.linalg.eigvalsh(a - b * np.eye(len(a)) - k / q)[0] > 0
    hi = min(float(data['projected_h'][i][i] / data['metric'][i][i]) for i in range(len(a)))
    lo, step = min(3 - eta, hi - 1), 1.
    for _ in range(64):
        if positive(lo):
            break
        lo -= step
        step *= 2
    else:
        raise ValueError('No lower proposal bracket')
    for _ in range(64):
        middle = (lo + hi) / 2
        if positive(middle):
            lo = middle
        else:
            hi = middle
    for check in range(8):
        margin = 10. ** (check - 9)
        lower = F(math.floor((lo - margin) * 10**10), 10**10)
        if pivots(data, lower) is not None:
            return lower, {'numerical_lower': lo, 'margin': margin, 'exact_checks': check + 1}
    raise ValueError('No exact lower bound accepted')


def witness_krylov_data(workspace, coefficients, steps=1):
    """Exact Rayleigh data on span{psi,...,H^steps psi}; upper only."""
    if type(steps) is not int or not 1 <= steps <= 3:
        raise ValueError('One to three upper-witness steps required')
    data = workspace_data(workspace)
    rayleigh(data, coefficients)  # validate the original physical witness
    model = workspace['model']
    psi = combine(*[(c, v) for c, v in zip(coefficients, workspace['basis'])])
    hpsi = combine(*[(c, v) for c, v in zip(coefficients, workspace['action'])])
    powers = [psi, hpsi]
    for _ in range(steps):
        v = powers[-1]
        powers.append(combine((1, model.reference_action(v)), (1, model.apply_polynomial(workspace['delta'], v))))
    vectors = powers[:-1]
    a = model.gram(vectors, powers[1:])
    if any(a[i][j] != a[j][i] for i in range(len(a)) for j in range(len(a))):
        raise ValueError('Nonsymmetric witness Hamiltonian')
    return {'metric': model.gram(vectors, vectors), 'projected_h': a, 'retained_dimension': len(vectors),
            'scope': 'Rayleigh upper-witness data only'}


def certificate_modes(certificate):
    modes = certificate.get('modes')
    if (certificate.get('kind') != 'fully_implicit_schur_v1' or type(modes) is not int
            or modes < 4 or modes % 2 or type(certificate.get('particles')) is not int
            or certificate['particles'] != modes // 2):
        raise ValueError('Unsupported fully implicit certificate')
    return modes


def targeting_guarantee(h, pairs, history, recipe=None, error_floor=F(1, 10**7), upper=None):
    from experiments.marginal_targeting_bound import context, ceiling, propose_retained_lower, check_step
    if not history:
        raise ValueError('Convergence certificate requires a target step')
    if recipe is None:
        initial = rayleigh(history[0][0], history[0][1])
        recipe = {'kind': 'residual_convergence_v1', 'upper': str(ceiling(initial, 10**9) if upper is None else F(upper)),
                  'error_floor': str(F(error_floor)), 'retained_lowers': []}
        generate = True
    else:
        if (type(recipe) is not dict or set(recipe) != {'kind', 'upper', 'error_floor', 'retained_lowers'}
                or recipe['kind'] != 'residual_convergence_v1'
                or type(recipe['retained_lowers']) is not list or len(recipe['retained_lowers']) != len(history)):
            raise ValueError('Malformed convergence recipe')
        generate = False
    guarantee = context(h, pairs, recipe['upper'], recipe['error_floor'])
    error = rayleigh(history[0][0], history[0][1]) - F(guarantee['ground_lower'])
    steps = []
    for i, (data, coefficients, next_data, next_coefficients) in enumerate(history):
        if generate:
            lower = propose_retained_lower(next_data, next_coefficients, guarantee)
            recipe['retained_lowers'].append(rational_text(lower))
        step = check_step(data, coefficients, next_data, next_coefficients, recipe['retained_lowers'][i], guarantee)
        error = F(guarantee['contraction_factor']) * error + F(step['additive_error_upper'])
        steps.append(step)
    return recipe, {'parameters': guarantee, 'steps': steps, 'iterated_energy_error_upper': rational_text(error),
                    'scope': 'Exact residual extensions and approximate Ritz error budgets verified; dimension and arithmetic cost are not bounded'}


def certificate_workspace(certificate, certify_floor=None):
    modes = certificate_modes(certificate)
    h = decode(certificate['hamiltonian'], modes, 4)
    workspace = enlarged_workspace(h, modes // 2)
    stages = [modes // 2 + 1, len(workspace['basis'])]
    target = certificate.get('target_coefficients')
    targets = certificate.get('target_chain')
    if targets is not None:
        if target is not None or type(targets) is not list or len(targets) > 3:
            raise ValueError('Expected an unambiguous chain of at most three targets')
    else:
        targets = [target] if target is not None else []
    history = []
    for i, coordinates in enumerate(targets):
        data = workspace_data(workspace)
        workspace = targeted_workspace(workspace, coordinates)
        stages.append(len(workspace['basis']))
        next_coefficients = targets[i + 1] if i + 1 < len(targets) else certificate.get('upper_coefficients')
        history.append((data, coordinates, workspace_data(workspace), next_coefficients))
    workspace['targeting_history'] = history
    recipe = certificate.get('targeting_guarantee')
    if recipe is not None or certify_floor is not None:
        if certificate.get('upper_krylov_coefficients') is not None:
            raise ValueError('Convergence requires the final vector in the retained workspace')
        if recipe is not None and certify_floor is not None:
            raise ValueError('Cannot replace an existing convergence recipe implicitly')
        recipe, receipt = targeting_guarantee(h, modes // 2, history, recipe,
                                            F(1, 10**7) if certify_floor is None else certify_floor)
        workspace['targeting_guarantee'] = recipe
        workspace['targeting_receipt'] = receipt
    return workspace, stages, targets


def certificate_upper(workspace, certificate):
    upper = rayleigh(workspace_data(workspace), certificate.get('upper_coefficients'))
    refinement = certificate.get('upper_krylov_coefficients')
    if refinement is not None:
        if type(refinement) is not list or not 2 <= len(refinement) <= 4:
            raise ValueError('Expected two to four upper-witness Krylov coefficients')
        upper = rayleigh(witness_krylov_data(workspace, certificate['upper_coefficients'], len(refinement) - 1), refinement)
    return upper


def replay(certificate):
    workspace, stages, _ = certificate_workspace(certificate)
    return workspace_receipt(certificate, workspace, stages)


def workspace_receipt(certificate, workspace, stages):
    modes = certificate_modes(certificate)
    data = workspace_data(workspace)
    lower = F(certificate['lower'])
    checked = pivots(data, lower)
    if checked is None:
        raise ValueError('Nonpositive exact lower certificate')
    upper = certificate_upper(workspace, certificate)
    if upper < lower:
        raise ValueError('Inconsistent interval')
    result = {'lower': rational_text(lower), 'upper': rational_text(upper), 'width': rational_text(upper - lower),
            'lower_float': float(lower), 'upper_float': float(upper), 'width_float': float(upper - lower),
            'modes': modes, 'particles': modes // 2,
            'stage_dimensions': stages,
            'basis_atom_count': data['basis_atom_count'], 'action_atom_count': data['action_atom_count'],
            'max_active_pairs': data['max_active_pairs'], 'schur_pivots': [rational_text(x) for x in checked],
            'scope': 'Both finite-model endpoints rebuilt in implicit coordinates; no general scaling claim'}
    if 'targeting_receipt' in workspace:
        result['targeting_convergence'] = workspace['targeting_receipt']
    return result


def certify_targeting(path, error_floor=F(1, 10**7)):
    source = Path(path)
    out = source.parent.with_name(source.parent.name + '_convergence')
    if out.exists():
        raise ValueError('Preserve previous convergence export')
    started = time.monotonic()
    certificate = json.loads(source.read_text())
    workspace, stages, _ = certificate_workspace(certificate, certify_floor=error_floor)
    receipt = workspace_receipt(certificate, workspace, stages)
    certificate['targeting_guarantee'] = workspace['targeting_guarantee']
    receipt['elapsed_seconds'] = time.monotonic() - started
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'path': str(out), 'stage_dimensions': stages, 'width_float': receipt['width_float'],
                      'error_floor': str(error_floor), 'elapsed_seconds': receipt['elapsed_seconds']}), flush=True)


def refine_upper(path, tag='', steps=1):
    if tag and (not tag.isidentifier() or not tag.isascii()):
        raise ValueError('Expected simple ASCII run tag')
    source = Path(path)
    if type(steps) is not int or not 1 <= steps <= 3:
        raise ValueError('One to three upper-witness steps required')
    out = source.parent.with_name(source.parent.name + '_upper' + str(steps) + ('_' + tag if tag else ''))
    if out.exists():
        raise ValueError('Preserve previous refinement')
    out.mkdir(parents=True)
    started = time.monotonic()
    certificate = json.loads(source.read_text())
    modes = certificate_modes(certificate)
    if certificate.get('upper_krylov_coefficients') is not None:
        raise ValueError('Expected unrefined fully implicit certificate')
    workspace, _, _ = certificate_workspace(certificate)
    data = workspace_data(workspace)
    lower = F(certificate['lower'])
    if pivots(data, lower) is None:
        raise ValueError('Original lower bound failed exact replay')
    original = rayleigh(data, certificate['upper_coefficients'])
    small = witness_krylov_data(workspace, certificate['upper_coefficients'], steps)
    coefficients, proposal = suggest_upper(small)
    upper = rayleigh(small, coefficients)
    if not lower <= upper <= original:
        raise ValueError('Refinement did not improve a consistent interval')
    certificate['upper_krylov_coefficients'] = coefficients
    # This recipe certifies a sequence of vectors in the retained spaces.
    # The refined energy certificate remains valid independently.
    certificate.pop('targeting_guarantee', None)
    (out / 'certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
    receipt = {'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
               'width_float': float(upper - lower), 'original_upper': str(original),
               'proposal': proposal, 'elapsed_seconds': time.monotonic() - started,
               'upper_steps': steps, 'scope': 'Upper-witness Hamiltonian steps; original lower certificate rechecked exactly'}
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('upper', 'width', 'original_upper', 'proposal')}), flush=True)


def enrich(path, error_floor=None):
    source = Path(path)
    out = source.parent.with_name(source.parent.name + '_enriched')
    if out.exists():
        raise ValueError('Preserve previous enrichment')
    out.mkdir(parents=True)
    started = time.monotonic()
    certificate = json.loads(source.read_text())
    workspace, stages, targets = certificate_workspace(certificate)
    if len(targets) >= 3:
        raise ValueError('Target-chain budget exhausted')
    data = workspace_data(workspace)
    old_lower = F(certificate['lower'])
    if pivots(data, old_lower) is None or certificate_upper(workspace, certificate) < old_lower:
        raise ValueError('Source interval failed exact validation')
    target, proposal = suggest_upper(data)
    history = workspace['targeting_history'][:]
    if history:
        history[-1] = (*history[-1][:3], target)
    previous_data = data
    (out / 'target_proposal.json').write_text(json.dumps(dict(proposal, coefficients=target), indent=2) + '\n')
    print(json.dumps({'stage': 'source_verified', 'dimensions': stages, 'seconds': time.monotonic() - started}), flush=True)
    workspace = targeted_workspace(workspace, target)
    stages.append(len(workspace['basis']))
    print(json.dumps({'stage': 'next_closure', 'dimensions': stages, 'seconds': time.monotonic() - started}), flush=True)
    data = workspace_data(workspace)
    coordinates = numerical_coordinates(data)
    upper_coordinates, proposal = suggest_upper(data, coordinates)
    lower, search = suggest_lower(data, coordinates)
    upper = rayleigh(data, upper_coordinates)
    if lower > upper:
        raise ValueError('Inconsistent interval')
    certificate.update(target_coefficients=None, target_chain=targets + [target],
                       upper_coefficients=upper_coordinates, lower=str(lower))
    certificate.pop('upper_krylov_coefficients', None)
    history.append((previous_data, target, data, upper_coordinates))
    previous_guarantee = certificate.get('targeting_guarantee', {})
    floor = error_floor if error_floor is not None else previous_guarantee.get('error_floor', F(1, 10**7))
    guarantee, convergence = targeting_guarantee(decode(certificate['hamiltonian'], certificate['modes'], 4),
                                                certificate['modes'] // 2, history, error_floor=floor,
                                                upper=previous_guarantee.get('upper'))
    certificate['targeting_guarantee'] = guarantee
    (out / 'certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
    receipt = {'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
               'width_float': float(upper - lower), 'stage_dimensions': stages,
               'proposal': proposal, 'search': search, 'targeting_convergence': convergence,
               'elapsed_seconds': time.monotonic() - started}
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('upper', 'width', 'proposal', 'search', 'targeting_convergence')}), flush=True)


def run(strength=F(1, 100), interaction=False, targeted=False, tag='', modes=10, error_floor=F(1, 10**7)):
    if type(modes) is not int or modes < 4 or modes % 2:
        raise ValueError('Even M>=4 required')
    if tag and (not tag.isidentifier() or not tag.isascii()):
        raise ValueError('Expected simple ASCII run tag')
    name = ('mixed' if interaction else 'cycle') + '_' + str(strength).replace('/', '_') + ('_targeted' if targeted else '')
    if modes != 10:
        name = 'm' + str(modes) + '_' + name
    if tag:
        name += '_' + tag
    out = Path(__file__).resolve().parents[1] / 'results/marginal_implicit_certificate' / name
    if out.exists():
        raise ValueError('Preserve previous run')
    out.mkdir(parents=True)
    started = time.monotonic()
    h = fixture(strength, interaction, modes)
    workspace = enlarged_workspace(h, modes // 2)
    data = workspace_data(workspace)
    coordinates = numerical_coordinates(data)
    upper_coefficients, proposal = suggest_upper(data, coordinates)
    (out / 'first_proposal.json').write_text(json.dumps(proposal, indent=2) + '\n')
    print(json.dumps({'stage': 'first_space', 'dimension': data['retained_dimension'], 'seconds': time.monotonic() - started,
                      **{k: v for k, v in proposal.items() if k not in ('upper', 'relative_squared_rounding_error')}}), flush=True)
    target = None
    if targeted:
        target = upper_coefficients
        first_data = data
        (out / 'target_coefficients.json').write_text(json.dumps(target) + '\n')
        workspace = targeted_workspace(workspace, target)
        print(json.dumps({'stage': 'target_closure', 'dimension': len(workspace['basis']), 'seconds': time.monotonic() - started}), flush=True)
        data = workspace_data(workspace)
        coordinates = numerical_coordinates(data)
        upper_coefficients, proposal = suggest_upper(data, coordinates)
    lower, search = suggest_lower(data, coordinates)
    upper = rayleigh(data, upper_coefficients)
    if upper < lower:
        raise ValueError('Inconsistent interval')
    certificate = {'kind': 'fully_implicit_schur_v1', 'modes': modes, 'particles': modes // 2,
                   'hamiltonian': encode(h), 'lower': str(lower), 'target_coefficients': target,
                   'upper_coefficients': upper_coefficients}
    convergence = None
    if targeted:
        guarantee, convergence = targeting_guarantee(h, modes // 2, [(first_data, target, data, upper_coefficients)],
                                                    error_floor=error_floor)
        certificate['targeting_guarantee'] = guarantee
    (out / 'certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
    receipt = {'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
               'width_float': float(upper - lower), 'retained_dimension': data['retained_dimension'],
               'elapsed_seconds': time.monotonic() - started, 'proposal': proposal, 'search': search,
               'scope': 'Exact accepted interval; independent replay required separately'}
    if convergence is not None:
        receipt['targeting_convergence'] = convergence
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('proposal', 'search', 'upper', 'width', 'targeting_convergence')}), flush=True)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify')
    p.add_argument('--refine-upper')
    p.add_argument('--enrich')
    p.add_argument('--certify-targeting')
    p.add_argument('--error-floor')
    p.add_argument('--upper-steps', type=int, choices=(1, 2, 3), default=1)
    p.add_argument('--strength', default='1/100')
    p.add_argument('--modes', type=int, default=10)
    p.add_argument('--interaction', action='store_true')
    p.add_argument('--targeted', action='store_true')
    p.add_argument('--tag', default='')
    args = p.parse_args()
    if args.certify_targeting:
        certify_targeting(args.certify_targeting, F(args.error_floor) if args.error_floor else F(1, 10**7))
    elif args.enrich:
        enrich(args.enrich, F(args.error_floor) if args.error_floor else None)
    elif args.refine_upper:
        refine_upper(args.refine_upper, args.tag, args.upper_steps)
    elif args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    else:
        run(F(args.strength), args.interaction, args.targeted, args.tag, args.modes,
            F(args.error_floor) if args.error_floor else F(1, 10**7))
