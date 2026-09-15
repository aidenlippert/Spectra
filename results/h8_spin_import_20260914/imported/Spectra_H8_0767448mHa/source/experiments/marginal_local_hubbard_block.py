"""Exact all-Fock reflection-block certificates for small Hubbard clusters.

Every particle/spin sector and both spatial-reflection characters are covered.
The supplied lower endpoint is checked against freshly constructed CAR actions.
"""
from fractions import Fraction as F
from math import lcm
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_polynomial_sos import integer_psd


def _exact(value, bound=10**6):
    if type(value) not in (int,str,F):
        raise ValueError('Bounded exact rational required')
    value = F(value)
    if abs(value) > bound or value.denominator > bound:
        raise ValueError('Rational magnitude or denominator exceeds bound')
    return value


def _check(L,U,t):
    if type(L) is not int or L not in (2,4,6):
        raise ValueError('L must be2,4,6')
    U,t = _exact(U),_exact(t)
    if U < 0 or t < 0:
        raise ValueError('Nonnegative Hubbard parameters required')
    return U,t


def _profiles(L,U,t,onsite=None,hopping=None):
    result = []
    for values,size,average in ((onsite,L,U),(hopping,L-1,t)):
        if values is None:
            values = [average]*size
        if type(values) is not list or len(values) != size:
            raise ValueError('Local profile has the wrong length')
        values = [_exact(x) for x in values]
        if any(x < 0 for x in values) or sum(values) != size*average:
            raise ValueError('Local profile must be nonnegative with its declared exact mean')
        if values != values[::-1]:
            raise ValueError('Local profile must preserve spatial reflection')
        result.append(values)
    return tuple(result)


def _density_profile(L,V=0,density=None):
    V=_exact(V)
    density=[V]*(L-1) if density is None else density
    if type(density) is not list or len(density)!=L-1:
        raise ValueError('Density profile has the wrong length')
    density=[_exact(x) for x in density]
    if sum(density)!=(L-1)*V or density!=density[::-1]:
        raise ValueError('Density profile must preserve reflection and its declared mean')
    return V,density


def _density_terms(L,density):
    terms=[]
    for i,V in enumerate(density):
        if not V:continue
        terms.append(((),V))
        for spin in (0,1):
            a,b=2*i+spin,2*(i+1)+spin
            terms.extend([(((1,a),(0,a)),-V),(((1,b),(0,b)),-V)])
            for other in (0,1):
                c=2*(i+1)+other
                terms.append((((1,a),(0,a),(1,c),(0,c)),V))
    return terms


def _terms(L,U,t,onsite=None,hopping=None,density=None):
    onsite = [U]*L if onsite is None else onsite
    hopping = [t]*(L-1) if hopping is None else hopping
    terms = []
    for i in range(L):
        terms.append((((1,2*i),(0,2*i),(1,2*i+1),(0,2*i+1)),onsite[i]))
    for i in range(L-1):
        for spin in (0,1):
            a,b = 2*i+spin,2*(i+1)+spin
            terms += [(((1,a),(0,b)),-hopping[i]),(((1,b),(0,a)),-hopping[i])]
    if density is not None:terms.extend(_density_terms(L,density))
    return [(word,a) for word,a in terms if a]


def _reflection(state,L):
    occupied = [2*(L-1-m//2)+m%2 for m in range(2*L) if state>>m&1]
    inversions = sum(a>b for j,a in enumerate(occupied) for b in occupied[j+1:])
    return sum(1<<m for m in occupied),(-1)**inversions


def _sector(state,L):
    up = sum((state>>(2*i))&1 for i in range(L))
    return up,state.bit_count()-up


def _actions(L,U,t,onsite=None,hopping=None,V=0,density=None):
    U,t = _check(L,U,t)
    onsite,hopping = _profiles(L,U,t,onsite,hopping)
    V,density = _density_profile(L,V,density)
    terms = _terms(L,U,t,onsite,hopping,density)
    actions = {}
    for state in range(4**L):
        shift = -sum(onsite[i]*(((state>>(2*i))&1)+((state>>(2*i+1))&1)-1)
                     for i in range(L))/2
        image = {state:shift}
        for word,a in terms:
            result = apply_word(word,state)
            if result:
                target,sign = result
                image[target] = image.get(target,0)+a*sign
        actions[state] = {s:a for s,a in image.items() if a}
    for state,image in actions.items():
        reflected,phase = _reflection(state,L)
        back,back_phase = _reflection(reflected,L)
        if back != state or phase*back_phase != 1:
            raise ValueError('Reflection must be an exact involution')
        reflected_image = {}
        for target,a in image.items():
            if _sector(target,L) != _sector(state,L):
                raise ValueError('Hamiltonian leaks from its spin-number sector')
            if actions[target].get(state,0) != a:
                raise ValueError('Exact CAR Hermiticity failed')
            r,s = _reflection(target,L)
            reflected_image[r] = a*s
        if reflected_image != {s:phase*a for s,a in actions[reflected].items()}:
            raise ValueError('Actual CAR action does not commute with reflection')
    return actions


def sector_matrices(sites,U,t,onsite=None,hopping=None,V=0,density=None):
    """Map (Nup,Ndown,parity) to (sector states, E^T K E, columns of E)."""
    U,t = _check(sites,U,t)
    actions = _actions(sites,U,t,onsite,hopping,V,density)
    sectors = {}
    for state in actions:
        sectors.setdefault(_sector(state,sites),[]).append(state)
    output = {}
    for key,states in sectors.items():
        seen = set()
        bases = {1:[],-1:[]}
        for state in states:
            if state in seen:
                continue
            reflected,phase = _reflection(state,sites)
            if reflected not in actions or _sector(reflected,sites) != key:
                raise ValueError('Reflection leaves its sector')
            seen.update((state,reflected))
            if reflected == state:
                bases[phase].append({state:1})
            else:
                for parity in (1,-1):
                    bases[parity].append({state:1,reflected:parity*phase})
        if sum(map(len,bases.values())) != len(states) or seen != set(states):
            raise ValueError('Incomplete reflection basis')
        for parity,columns in bases.items():
            if not columns:
                continue
            lookup = {}
            for i,column in enumerate(columns):
                norm = sum(a*a for a in column.values())
                if norm not in (1,2):
                    raise ValueError('Invalid orbit Gram norm')
                for state,a in column.items():
                    if state in lookup:
                        raise ValueError('Reflection columns overlap within a character')
                    lookup[state] = (i,a)
            matrix = [[F(0) for _ in columns] for _ in columns]
            for j,column in enumerate(columns):
                for source,a in column.items():
                    for target,b in actions[source].items():
                        if target in lookup:
                            i,c = lookup[target]
                            matrix[i][j] += a*b*c
            if any(matrix[i][j] != matrix[j][i] for i in range(len(matrix)) for j in range(i)):
                raise ValueError('Projected Hermiticity failed')
            output[(*key,parity)] = (states,matrix,columns)
    if sum(len(item[1]) for item in output.values()) != 4**sites:
        raise ValueError('All-Fock dimension coverage failed')
    return output


def replay(certificate):
    if type(certificate) is not dict or certificate.get('kind') != 'local_hubbard_block_v1':
        raise ValueError('Unsupported local certificate')
    L = certificate.get('sites')
    U,t = _check(L,certificate.get('U'),certificate.get('t'))
    onsite,hopping = _profiles(L,U,t,certificate.get('onsite_profile'),certificate.get('hopping_profile'))
    V,density = _density_profile(L,certificate.get('V',0),certificate.get('density_profile'))
    lower = _exact(certificate.get('lower'),10**9)
    data = sector_matrices(L,U,t,onsite,hopping,V,density)
    sectors = []
    for key,(_,matrix,columns) in data.items():
        dimension = len(matrix)
        norms = [sum(a*a for a in column.values()) for column in columns]
        shifted = [[matrix[i][j]-(lower*norms[i] if i==j else 0)
                    for j in range(dimension)] for i in range(dimension)]
        scale = lcm(*(x.denominator for row in shifted for x in row))
        integer = [[int(x*scale) for x in row] for row in shifted]
        result = integer_psd(integer)
        sectors.append({'sector':list(key),'dimension':dimension,'psd':result})
    norm = quotient = None
    if 'upper_vector' in certificate:
        submitted = certificate['upper_vector']
        if type(submitted) is not dict or not 1 <= len(submitted) <= 4**L:
            raise ValueError('Bounded nonempty physical upper vector required')
        vector = {}
        for key,a in submitted.items():
            if type(key) not in (int,str):
                raise ValueError('Integer determinant label required')
            state = int(key)
            if (state in vector or not 0 <= state < 4**L or state.bit_count() != L
                    or type(a) is not int or abs(a) > 10**12):
                raise ValueError('Invalid half-filled exact upper vector')
            vector[state] = a
        norm = sum(a*a for a in vector.values())
        if norm <= 0:
            raise ValueError('Upper vector must have positive norm')
        # For nonuniform profiles the centered shift need not vanish even at
        # half filling. Its actual expectation must be retained.
        energy = F(0)
        for source,a in vector.items():
            shift = -sum(onsite[i]*(((source>>(2*i))&1)+((source>>(2*i+1))&1)-1)
                         for i in range(L))/2
            energy += a*a*shift
            for word,b in _terms(L,U,t,onsite,hopping,density):
                image = apply_word(word,source)
                if image:
                    energy += a*b*image[1]*vector.get(image[0],0)
        quotient = energy/norm
        if quotient < lower:
            raise ValueError('Local lower exceeds physical upper')
    return {'accepted':True,'sites':L,'U':str(U),'t':str(t),'V':str(V),'sectors':sectors,
            'onsite_profile':list(map(str,onsite)),'hopping_profile':list(map(str,hopping)),
            'density_profile':list(map(str,density)),
            'sector_count':len(sectors),'sum_dimensions':sum(x['dimension'] for x in sectors),
            'maximum_psd_dimension':max(x['dimension'] for x in sectors),
            'lower':str(lower),'upper_norm':str(norm) if norm is not None else None,
            'upper_quotient':str(quotient) if quotient is not None else None,
            'reflection_commutes_with_original_car':True,
            'scope':'Exact centered local Hubbard operator with optional nearest-neighbor density interaction, positive in every Fock sector after subtracting the lower endpoint. Dense reflected matrices remain bounded to clusters of2,4,6sites; no global sector assumption.'}
