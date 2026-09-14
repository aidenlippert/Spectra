"""Exact dual limits for a fixed joint-projector and diagonal-correction span."""
from fractions import Fraction as F

from experiments.marginal_joint_family_limit import _physical_vector, _weight, _profiles
from experiments.marginal_local_hubbard_block import _actions, _exact, _sector
from experiments.marginal_projector_extendibility import _projector_vector, _telescoping_diagonal
from experiments.marginal_charged_projectors import charged_vectors


def replay(c):
    if type(c) is not dict or c.get('kind') != 'joint_diagonal_family_limit_v1':
        raise ValueError('Unsupported diagonal family limit')
    shapes = c.get('diagonal_shapes')
    if type(shapes) is not list or not 1 <= len(shapes) <= 9:
        raise ValueError('One through nine diagonal shapes required')
    shapes = [_telescoping_diagonal(shape) for shape in shapes]
    if len(set().union(*(set(shape) for shape in shapes))) > 64:
        raise ValueError('Diagonal span exceeds the 64-entry energy certificate cap')
    half, hn = _projector_vector(c.get('half_vector'), 6)
    if {_sector(s, 6) for s, a in half.items() if a} != {(3, 3)}:
        raise ValueError('Half projector must occupy spin sector (3,3)')
    charged, cn = charged_vectors(c.get('charged_vector'))
    ratio = _exact(c.get('ratio'), 10**9)
    th, tj = (_exact(c.get(name), 10**12) for name in ('theta_half', 'theta_joint'))
    if ratio <= 0 or not 0 <= th <= 1 or not 0 <= tj <= max(1, ratio):
        raise ValueError('Invalid fixed-projector ceiling or ratio')
    mixture = c.get('mixture')
    # Trace, six profile moments, k diagonal moments and two fidelity rows
    # give 9+k scalar constraints in the finite dual LP (including slacks).
    if type(mixture) is not list or not 1 <= len(mixture) <= 9+len(shapes):
        raise ValueError('Mixture exceeds the bounded scalar constraint count')
    states = []
    for item in mixture:
        if type(item) is not dict:
            raise ValueError('Explicit physical mixture required')
        vector, norm = _physical_vector(item.get('vector'))
        states.append((_weight(item.get('weight')), vector, norm))
    if sum(w for w, _, _ in states) != 1:
        raise ValueError('Mixture trace must equal one exactly')
    moments = [sum((w*sum((a*a*(shape.get(s & 1023, F(0))-shape.get(s >> 2, F(0)))
                          for s, a in vector.items()), F(0))/norm
                    for w, vector, norm in states), F(0)) for shape in shapes]
    if any(moments):
        raise ValueError('Every added diagonal expectation must cancel exactly')
    def fidelity(source, source_norm):
        return sum((w*F(sum(a*source.get(s, 0) for s, a in vector.items())**2,
                        norm*source_norm) for w, vector, norm in states), F(0))
    ph = fidelity(half, hn)
    q = ph+ratio*sum((fidelity(v, cn) for v in charged), F(0))
    if ph > th or q > tj:
        raise ValueError('Mixture violates a fixed-projector fidelity ceiling')
    energies = []
    for index in range(7):
        x = [F(int(index == i+1)) for i in range(6)]
        onsite, hopping, density = _profiles(x)
        actions = _actions(6, F(10, 3), 1, onsite, hopping, F(1, 2), density)
        energies.append(sum((w*sum((F(a)*b*vector.get(target, 0)
                                    for s, a in vector.items() for target, b in actions[s].items()), F(0))/norm
                             for w, vector, norm in states), F(0)))
    gradients = [energy-energies[0] for energy in energies[1:]]
    if any(gradients):
        raise ValueError('Profile expectations must cancel exactly')
    upper = energies[0]/5
    if 'proposed_periodic_family_upper' in c and F(c['proposed_periodic_family_upper']) != upper:
        raise ValueError('Proposed family limit disagrees with physical expectation')
    return {'accepted': True, 'periodic_family_upper': str(upper),
            'mixture_sources': len(states), 'mixture_trace': '1',
            'profile_gradients': list(map(str, gradients)), 'diagonal_moments': list(map(str, moments)),
            'diagonal_shapes': len(shapes), 'diagonal_union_entries': len(set().union(*(set(s) for s in shapes))),
            'half_expectation': str(ph), 'joint_expectation': str(q),
            'theta_half': str(th), 'theta_joint': str(tj),
            'scope': 'Upper limit on attainable periodic LOWER certificates for fixed six-site U4,t1,V1/2 projector sources, ratio, ceilings and the supplied diagonal span. All reflected mean-correct profiles and unrestricted real coefficients of these shapes are covered. Nonnegative projector penalties are required. This local PSD mixture need not extend globally and is not a physical ground-energy upper bound. Projector ceiling validity is certified separately. Other diagonal shapes, changed sources or ratio, and general representability are not covered.'}
