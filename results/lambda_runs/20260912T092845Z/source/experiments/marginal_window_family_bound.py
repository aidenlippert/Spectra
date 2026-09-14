"""Exact primal/dual bracket for a two-parameter four-site window family.

For every normalized physical rho, lambda_min(K(a,b)) <= tr(rho K(a,b)).
The right side is affine; its maximum over the declared rectangle occurs
at a corner. A nearly stationary rational Rayleigh witness bounds the entire
family without relying on the numerical optimizer or a PSD relaxation dual.
"""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _actions, _exact, replay as replay_local


def local_certificate(a,b,lower):
    a,b = _exact(a),_exact(b)
    if not 0 <= a <= 6 or not 0 <= b <= F(3,2):
        raise ValueError('Window parameters outside the declared rectangle')
    return {'kind':'local_hubbard_block_v1','sites':4,'U':'3','t':'1',
            'onsite_profile':list(map(str,[a,6-a,6-a,a])),
            'hopping_profile':list(map(str,[b,3-2*b,b])), 'lower':str(lower)}


def _vector(vector):
    if type(vector) is not dict or not 1 <= len(vector) <= 256:
        raise ValueError('Bounded nonempty physical vector required')
    amplitudes = {}
    for key,value in vector.items():
        if type(key) not in (int,str):
            raise ValueError('Integer state label required')
        state = int(key)
        if (state in amplitudes or not 0 <= state < 256 or state.bit_count() != 4
                or type(value) is not int or abs(value) > 10**12):
            raise ValueError('Bounded half-filled integer vector required')
        amplitudes[state] = value
    norm = sum(x*x for x in amplitudes.values())
    if norm <= 0:
        raise ValueError('Positive vector norm required')
    return amplitudes,norm


def affine_rayleigh(vector):
    amplitudes,norm = _vector(vector)
    energies = []
    for a,b in ((0,0),(1,0),(0,1)):
        actions = _actions(4,3,1,[a,6-a,6-a,a],[b,3-2*b,b])
        energies.append(sum(value*coefficient*amplitudes.get(target,0)
                            for source,value in amplitudes.items()
                            for target,coefficient in actions[source].items())/F(norm))
    return norm,(energies[0],energies[1]-energies[0],energies[2]-energies[0])


def overlap_consistency(vector):
    """Exact three-site reductions of the physical four-site pure witness."""
    amplitudes,norm = _vector(vector)
    reductions = []
    for right in (False,True):
        groups = {}
        for state,a in amplitudes.items():
            environment = state&3 if right else state>>6
            local = state>>2 if right else state&63
            groups.setdefault(environment,{})[local] = a
        matrix = {}
        for part in groups.values():
            for i,a in part.items():
                for j,b in part.items():
                    matrix[i,j] = matrix.get((i,j),0)+a*b
        reductions.append({k:x for k,x in matrix.items() if x})
    left,right = reductions
    equal = left == right
    purity = F(sum(x*x for x in left.values()),norm*norm)
    if not 0 < purity <= 1:
        raise ValueError('Invalid reduced-state purity')
    return {'three_site_reductions_equal':equal,'left_nonzero_entries':len(left),
            'right_nonzero_entries':len(right),'three_site_purity':str(purity),
            'no_five_site_translation_invariant_extension':equal and purity < 1,
            'extension_obstruction':'A five-site extension of a pure four-site marginal must factor across that four-site block and the last site. Its shifted four-site marginal then has purity at most the three-site purity; if this is below1 it cannot equal the original pure marginal.'}


def replay(certificate):
    if type(certificate) is not dict or certificate.get('kind') != 'hubbard_four_window_family_bound_v1':
        raise ValueError('Unsupported window-family certificate')
    local = local_certificate(certificate.get('a'),certificate.get('b'),certificate.get('lower'))
    local_receipt = replay_local(local)
    norm,(constant,alpha,beta) = affine_rayleigh(certificate.get('upper_vector'))
    consistency = overlap_consistency(certificate.get('upper_vector'))
    if consistency['three_site_reductions_equal'] and (alpha or beta):
        raise ValueError('Translation consistency disagrees with affine coefficients')
    corners = [constant+alpha*a+beta*b for a in (0,6) for b in (0,F(3,2))]
    upper = max(corners)
    lower = F(local_receipt['lower'])
    if lower > upper:
        raise ValueError('Local lower exceeds family-wide upper')
    if 'target_family_upper' in certificate:
        if upper > _exact(certificate['target_family_upper'],10**12):
            raise ValueError('Affine witness does not prove the requested family upper')
    return {'accepted':True,'local_minimum_lower':str(lower),
            'maximum_local_minimum_upper':str(upper),'family_bracket_width':str(upper-lower),
            'norm':str(norm),'affine_expectation':list(map(str,(constant,alpha,beta))),
            'corner_expectations':list(map(str,corners)),'local_certificate':local,
            'local_replay':local_receipt,
            'overlap_consistency':consistency,
            'arbitrary_three_site_boundary_correction_upper':str(constant) if consistency['three_site_reductions_equal'] else None,
            'boundary_correction_scope':'When the two three-site reductions agree, every Hermitian even correction A_left-A_right has zero expectation. The constant upper then applies to the uniform centered four-site window plus any such correction, not only to the two-parameter profile family.',
            'scope':'Exact bound on max_{a in[0,6],b in[0,3/2]} lambda_min of the centered four-site window with onsite[a,6-a,6-a,a] and hopping[b,3-2b,b]. This upper limits the strength of that lower-certificate family; it is not an upper bound on the full-chain ground energy.'}
