"""Exact dual limits for a fixed joint-projector and diagonal-correction span."""
from fractions import Fraction as F

from experiments.marginal_joint_family_limit import _physical_vector, _weight, _profiles
from experiments.marginal_local_hubbard_block import _actions, _exact, _sector
from experiments.marginal_projector_extendibility import _projector_vector, _telescoping_diagonal
from experiments.marginal_charged_projectors import charged_vectors


def replay(c):
    return _replay(c, free_range_two_profile=False)


def _replay(c, *, free_range_two_profile, full_quadratic_charge=False, charge_square_pairs=False, full_charge_indicators=False, full_signed_charge=False, hopping_telescope=False, spin_telescope=False, spectator_hopping=False):
    """Shared exact constraints with bounded, physically specified extra rows."""
    if type(c) is not dict or c.get('kind') != 'joint_diagonal_family_limit_v1':
        raise ValueError('Unsupported diagonal family limit')
    if 'W' in c or 'range_two_density_profile' in c:
        raise ValueError('Range-two fields require the explicit range-two family verifier')
    if 'quadratic_charge_telescope' in c or 'charge_square_pair_telescope' in c or 'higher_charge_indicator_telescope' in c or 'signed_charge_telescope' in c or 'hopping_telescope' in c or 'spin_telescope' in c or 'spectator_hopping' in c or 'pair_transfer' in c:
        raise ValueError('Family limits constrain moments, not fixed quadratic coefficients')
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
    if type(free_range_two_profile) is not bool or type(full_quadratic_charge) is not bool or (full_quadratic_charge and not free_range_two_profile):
        raise ValueError('Explicit bounded profile mode required')
    if type(charge_square_pairs) is not bool or (charge_square_pairs and not full_quadratic_charge):
        raise ValueError('Charge-square family requires full quadratic mode')
    if type(full_charge_indicators) is not bool or (full_charge_indicators and not charge_square_pairs):
        raise ValueError('Full indicator family requires charge-square pairs')
    if type(full_signed_charge) is not bool or (full_signed_charge and not full_charge_indicators):
        raise ValueError('Full signed-charge family requires full indicator mode')
    if type(hopping_telescope) is not bool or (hopping_telescope and not full_signed_charge):
        raise ValueError('Hopping family requires full signed-charge mode')
    if type(spin_telescope) is not bool or (spin_telescope and not hopping_telescope):
        raise ValueError('Spin family requires hopping telescope mode')
    if type(spectator_hopping) is not bool or (spectator_hopping and not spin_telescope):
        raise ValueError('Spectator hopping family requires spin telescope mode')
    if type(mixture) is not list or not 1 <= len(mixture) <= 9+len(shapes)+int(free_range_two_profile)+int(full_quadratic_charge)+4*int(charge_square_pairs)+6*int(full_charge_indicators)+36*int(full_signed_charge)+int(hopping_telescope)+4*int(spin_telescope)+14*int(spectator_hopping):
        raise ValueError('Mixture exceeds the bounded scalar constraint count')
    states = []
    for item in mixture:
        if type(item) is not dict:
            raise ValueError('Explicit physical mixture required')
        vector, norm = _physical_vector(item.get('vector'))
        states.append((_weight(item.get('weight')), vector, norm))
    if sum(w for w, _, _ in states) != 1:
        raise ValueError('Mixture trace must equal one exactly')
    range_two_gradient = F(0)
    if free_range_two_profile:
        from experiments.marginal_range_two_density import diagonal_value
        range_two_gradient = sum((w*sum((a*a*diagonal_value(s, [1, -1, -1, 1])
                                        for s, a in vector.items()), F(0))/norm
                                  for w, vector, norm in states), F(0))
        if range_two_gradient:
            raise ValueError('Free range-two profile expectation must cancel exactly')
    quadratic_moments = {}
    if full_quadratic_charge:
        from experiments.marginal_quadratic_charge_telescope import LABELS, local_value
        # Five of these six directions are implied by the profile constraints.
        # Check all six exactly; the remaining independent row adds one source.
        for label in LABELS:
            moment = sum((w*sum((a*a*local_value(s, {label: F(1)})
                                for s, a in vector.items()), F(0))/norm
                          for w, vector, norm in states), F(0))
            if moment:
                raise ValueError('Every quadratic charge expectation must cancel exactly')
            quadratic_moments[label] = str(moment)
    square_moments = {}
    if charge_square_pairs:
        from experiments.marginal_charge_square_pairs import LABELS, local_value
        for label in LABELS:
            moment = sum((w*sum((a*a*local_value(s, {label: F(1)})
                                for s, a in vector.items()), F(0))/norm
                          for w, vector, norm in states), F(0))
            if moment:
                raise ValueError('Every charge-square pair expectation must cancel exactly')
            square_moments[label] = str(moment)
    indicator_moments = {}
    if full_charge_indicators:
        from experiments.marginal_charge_indicator_telescope import ALL_LABELS, local_value
        # Check the entire indicator space, including the six old directions.
        for label in ALL_LABELS:
            moment = sum((w*sum((a*a*local_value(s, {label: F(1)})
                                for s, a in vector.items()), F(0))/norm
                          for w, vector, norm in states), F(0))
            if moment:
                raise ValueError('Every charge indicator expectation must cancel exactly')
            indicator_moments[label] = str(moment)
    signed_moments = {}
    if full_signed_charge:
        from experiments.marginal_signed_charge_telescope import PATTERNS, component
        signed_moments = dict.fromkeys(PATTERNS, F(0))
        for w,vector,norm in states:
            numerators = dict.fromkeys(PATTERNS, 0)
            for s,a in vector.items():
                for five,orientation in [(s&1023,1),(s>>2,-1)]:
                    entry = component(five)
                    if entry is not None:
                        label,sign = entry
                        numerators[label] += orientation*sign*a*a
            for label,numerator in numerators.items():
                signed_moments[label] += w*F(numerator,norm)
        if any(signed_moments.values()):
            raise ValueError('Every signed-charge pattern expectation must cancel exactly')
    hopping_moment = F(0)
    if hopping_telescope:
        from experiments.marginal_hopping_telescope import actions
        action = actions()
        hopping_moment = sum((w*F(sum(a*b*vector.get(t,0) for s,a in vector.items() for t,b in action[s].items()),norm) for w,vector,norm in states),F(0))
        if hopping_moment:
            raise ValueError('Hopping telescope expectation must cancel exactly')
    spin_moments = {}
    if spin_telescope:
        from experiments.marginal_spin_telescope import LABELS, actions
        for label in LABELS:
            action = actions({label:1})
            moment = sum((w*F(sum(a*b*vector.get(t,0) for s,a in vector.items() for t,b in action[s].items()),norm) for w,vector,norm in states),F(0))
            if moment:raise ValueError('Every spin telescope expectation must cancel exactly')
            spin_moments[label] = str(moment)
    spectator_moments = {}
    if spectator_hopping:
        from experiments.marginal_spectator_hopping import LABELS, actions
        for label in LABELS:
            action = actions({label:1})
            moment = sum((w*F(sum(a*b*vector.get(t,0) for s,a in vector.items() for t,b in action[s].items()),norm) for w,vector,norm in states),F(0))
            if moment:raise ValueError('Every spectator hopping expectation must cancel exactly')
            spectator_moments[label] = str(moment)
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
    result = {'accepted': True, 'periodic_family_upper': str(upper),
            'mixture_sources': len(states), 'mixture_trace': '1',
            'profile_gradients': list(map(str, gradients)), 'diagonal_moments': list(map(str, moments)),
            'diagonal_shapes': len(shapes), 'diagonal_union_entries': len(set().union(*(set(s) for s in shapes))),
            'half_expectation': str(ph), 'joint_expectation': str(q),
            'theta_half': str(th), 'theta_joint': str(tj),
            'scope': 'Upper limit on attainable periodic LOWER certificates for fixed six-site U4,t1,V1/2 projector sources, ratio, ceilings and the supplied diagonal span. All reflected mean-correct profiles and unrestricted real coefficients of these shapes are covered. Nonnegative projector penalties are required. This local PSD mixture need not extend globally and is not a physical ground-energy upper bound. Projector ceiling validity is certified separately. Other diagonal shapes, changed sources or ratio, and general representability are not covered.'}
    if free_range_two_profile:
        result['range_two_profile_gradient'] = str(range_two_gradient)
    if full_quadratic_charge:
        result['quadratic_charge_moments'] = quadratic_moments
    if charge_square_pairs:
        result['charge_square_pair_moments'] = square_moments
    if full_charge_indicators:
        result['charge_indicator_moments'] = indicator_moments
    if full_signed_charge:
        result['signed_charge_moments'] = {key:str(value) for key,value in signed_moments.items()}
    if hopping_telescope:
        result['hopping_telescope_moment'] = str(hopping_moment)
    if spin_telescope:
        result['spin_telescope_moments'] = spin_moments
    if spectator_hopping:
        result['spectator_hopping_moments'] = spectator_moments
    return result
