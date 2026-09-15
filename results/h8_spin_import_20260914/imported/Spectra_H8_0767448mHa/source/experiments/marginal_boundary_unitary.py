"""Bounded endpoint-RDM contractions for a physical contact-unitary circuit."""
from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce

from experiments.marginal_local_hubbard_block import _exact, _actions
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
from experiments.marginal_tiled_upper import _uniform_h, replay_upper8, remainder_upper
from experiments.marginal_polynomial_sos import integer_psd
from results.marginal_graded_hubbard8.discovery.singlet_moment_energy import upper_monomials


def _physical_source(c,h):
    """Validate and reconstruct the shared original H8 orbit-state recipe."""
    checked = replay_upper8(c,h)
    oracle = SymmetryMomentOracle(_uniform_h(8))
    vectors = [oracle.compress({int(s):a for s,a in v.items()}) for v in c['embedding']['basis']]
    state = {}
    for row in reversed(upper_monomials(c['upper_chebyshev_coefficients'])):
        state = oracle.apply(state)
        for a,v in zip(row,vectors):
            for s,b in v.items():state[s] = state.get(s,0)+a*b
        state = {s:a for s,a in state.items() if a}
        if len(state)>4096:raise ValueError('Existing orbit-support budget exceeded')
    if not state:raise ValueError('Nonzero physical trial state required')
    denominator = lcm(*(a.denominator for a in state.values()))
    state = {s:int(a*denominator) for s,a in state.items()}
    divisor = reduce(gcd,(abs(a) for a in state.values()))
    state = {s:a//divisor for s,a in state.items()}
    norm = oracle.dot(state,state)
    energy = F(oracle.dot(state,oracle.apply(state)),norm)
    if energy != F(checked['upper']):raise ValueError('Independent physical energy disagreement')
    return oracle,state,checked


def _state_reductions(oracle,state,checked=None):
    """Read a physical orbit state constructed by one of the bounded callers."""
    if not state or len(state)>4096:raise ValueError('Bounded nonzero orbit state required')
    for r in state:
        entry=oracle.orbit(r)
        if entry is None or entry[0]!=r:raise ValueError('Canonical physical orbit required')
    norm=oracle.dot(state,state)
    if norm<=0:raise ValueError('Positive physical state norm required')
    energy=F(oracle.dot(state,oracle.apply(state)),norm)
    total_doublons,total_density = 0,0
    for r,a in state.items():
        for s,phase in oracle.orbit(r)[3].items():
            weight=a*a
            charges=[((s>>(2*i))&3).bit_count()-1 for i in range(8)]
            total_doublons+=weight*sum(((s>>(2*i))&3)==3 for i in range(8))
            total_density+=weight*sum(charges[i]*charges[i+1] for i in range(7))
    matrices = []
    for offset in (0,12):
        matrix = [[0]*16 for _ in range(16)]
        mask = 15<<offset
        for r,a in state.items():
            for s,phase in oracle.orbit(r)[3].items():
                row = (s>>offset)&15
                for col in range(16):
                    target = (s&~mask)|(col<<offset)
                    if not oracle.valid_state(target):continue
                    entry = oracle.orbit(target)
                    if entry is not None:
                        matrix[row][col] += a*phase*state.get(entry[0],0)*entry[1]
        if sum(matrix[i][i] for i in range(16)) != norm:
            raise ValueError('Edge density trace must equal the physical norm')
        if any(matrix[i][j]!=matrix[j][i] for i in range(16) for j in range(16)):
            raise ValueError('Edge density Hermiticity failed')
        matrices.append([[F(a,norm) for a in row] for row in matrix])
    return {'energy8':str(energy),'norm':str(norm),'orbit_support':len(state),
            'total_doublons':str(F(total_doublons,norm)),
            'internal_density_correlation':str(F(total_density,norm)),
            'left':[[str(a) for a in row] for row in matrices[0]],
            'right':[[str(a) for a in row] for row in matrices[1]],
            'upper_replay':checked,
            'scope':'Exact physical two-site partial traces, streamed from a bounded caller-constructed H8 orbit state. No supplied RDM is accepted as physical.'}


def edge_reductions(c,h):
    """Recompute two physical16x16 edge reductions from the bounded H8 state."""
    oracle,state,checked=_physical_source(c,h)
    return _state_reductions(oracle,state,checked)


def rotation(x):
    """Real rational unitary on normalized covalent/ionic two-site singlets."""
    x = _exact(x,10**6)
    cosine,sine = (1-x*x)/(1+x*x),2*x/(1+x*x)
    singlet,ionic = {9:1,6:-1},{3:1,12:1}
    columns = {}
    for j in range(16):
        col = {}
        for i in range(16):
            p = F(singlet.get(i,0)*singlet.get(j,0)+ionic.get(i,0)*ionic.get(j,0),2)
            a = F(ionic.get(i,0)*singlet.get(j,0)-singlet.get(i,0)*ionic.get(j,0),2)
            value = F(i==j)+(cosine-1)*p+sine*a
            if value:col[i]=value
        columns[j]=col
    for i in range(16):
        for j in range(16):
            if sum(a*columns[j].get(s,0) for s,a in columns[i].items()) != int(i==j):
                raise ValueError('Exact contact unitarity failed')
        spin = ((i&5).bit_count(),(i&10).bit_count())
        if any(((s&5).bit_count(),(s&10).bit_count())!=spin for s in columns[i]):
            raise ValueError('Contact rotation must preserve both spin populations')
    return columns


def _embedded_gate(columns):
    return {s:{(s&195)|(t<<2):a for t,a in columns[(s>>2)&15].items()}
            for s in range(256)}


def _density(left,right):
    matrices = []
    for matrix in (left,right):
        if type(matrix) is not list or len(matrix)!=16 or any(type(row) is not list or len(row)!=16 for row in matrix):
            raise ValueError('Two16x16 edge matrices required')
        if any(type(a) not in (int,str,F) for row in matrix for a in row):
            raise ValueError('Exact edge matrix entries required')
        parsed = [[F(a) for a in row] for row in matrix]
        if any(parsed[i][j]!=parsed[j][i] for i in range(16) for j in range(16)):
            raise ValueError('Real symmetric edge matrices required')
        matrices.append(parsed)
    a,b = matrices
    return lambda i,j:a[i&15][j&15]*b[i>>4][j>>4]


def _operator_shift(data,columns,norm_factor,U=4,t=1,V=0):
    density = _density(data['right'],data['left'])
    gate = _embedded_gate(columns)
    # Restore the physical operator from the centered local engine. For a
    # nonunitary filter, commutation with N alone would not justify dropping
    # a centered shift after normalization. All256 local states are included.
    actions = _actions(4,U,t)
    for s in actions:
        charges=[((s>>(2*i))&3).bit_count()-1 for i in range(4)]
        actions[s][s]=actions[s].get(s,0)+F(U,2)*(s.bit_count()-4)+V*sum(charges[i]*charges[i+1] for i in range(3))
    before = sum(a*density(s,t) for s,image in actions.items() for t,a in image.items())
    transformed = {}
    for s,column in gate.items():
        image = {}
        for k,a in column.items():
            for t,b in actions[k].items():image[t]=image.get(t,0)+a*b
        transformed[s]=image
    after = F(0)
    for i,column in gate.items():
        for j,image in transformed.items():
            rho = density(i,j)
            if rho:
                after += rho*sum(a*image.get(s,0) for s,a in column.items())
    return after/norm_factor-before


def energy_shift(data,x):
    """Algebraic contraction; replay always supplies freshly physical data."""
    return _operator_shift(data,rotation(x),F(1))


def model_shift(data,parameter,operation='unitary',U=4,t=1,V=0):
    """Boundary shift for a separately specified nearest-neighbor target."""
    U,t,V = _exact(U),_exact(t),_exact(V)
    if U<0 or t<0:raise ValueError('Nonnegative target U,t required')
    parameter = _exact(parameter)
    if operation=='unitary':
        return _operator_shift(data,rotation(parameter),F(1),U,t,V)
    if operation!='linear_filter':raise ValueError('Unsupported boundary operation')
    # The filter's endpoint-preserving channel identity requires exact half
    # occupation of both spin modes on both original contacts.
    for name,site in (('right',1),('left',0)):
        diagonal=[F(data[name][i][i]) for i in range(16)]
        for spin in (0,1):
            if sum(a*((i>>(2*site+spin))&1) for i,a in enumerate(diagonal))!=F(1,2):
                raise ValueError('Linear filter requires half-occupied contact spins')
    h = _actions(2,F(0),F(1))
    columns={i:{j:F(i==j)-parameter*h[i].get(j,0) for j in range(16)
                if F(i==j)-parameter*h[i].get(j,0)} for i in range(16)}
    return _operator_shift(data,columns,1+parameter*parameter,U,t,V)


def replay_model(c,h,parameter,sites,operation='unitary',U=4,t=1,V=0):
    """Physical source state is fixed; target U,t,V energies are recomputed."""
    if type(sites) is not int or sites%8 or not 8<=sites<=10**9:
        raise ValueError('Transferred target currently requires multiples of8 sites')
    U,t,V = _exact(U),_exact(t),_exact(V)
    if U<0 or t<0:raise ValueError('Nonnegative target U,t required')
    parameter = _exact(parameter)
    if operation not in ('unitary','linear_filter'):raise ValueError('Unsupported boundary operation')
    data=edge_reductions(c,h)
    doublons=F(data['total_doublons'])
    # Source H8=4D-T is validated independently, so its full kinetic
    # expectation follows exactly without another long-range contraction.
    kinetic=4*doublons-F(data['energy8'])
    block=U*doublons-t*kinetic+V*F(data['internal_density_correlation'])
    charges=[]
    for name,site in (('right',1),('left',0)):
        charges.append(sum(F(data[name][i][i])*(((i>>(2*site))&3).bit_count()-1)
                           for i in range(16)))
    baseline=V*charges[0]*charges[1]
    shift=model_shift(data,parameter,operation,U,t,V)
    q=sites//8
    upper=q*block+(q-1)*(baseline+shift)
    return {'accepted':True,'sites':sites,'blocks':q,'operation':operation,
            'parameter':str(parameter),'target':{'U':str(U),'t':str(t),'V':str(V)},
            'target_block_energy':str(block),'unfiltered_cut_energy':str(baseline),
            'gate_energy_shift':str(shift),'upper':str(upper),'upper_per_site':str(upper/sites),
            'edge_reductions':data,
            'scope':'Exact physical upper for H=-t*sum hopping+U*sum doublons+V*sum(n_i-1)(n_{i+1}-1). The source trial state is generated from the validated original H8 U4 t1 operator, and its target energy is freshly recomputed. Contact operations preserve remote endpoint RDMs, so shifts add. No target ground-energy lower or target-state optimality claim.'}


def shift_polynomial(data):
    """Exact numerator of the shift over (1+x²)²; degree is at most four.

    The gate numerator is I+2x*A+x²*(I-2P), so multiplying its two
    occurrences around H proves the degree bound before interpolation.
    """
    points = (0,1,-1,2,-2)
    rows = [[F(x)**k for k in range(5)]+[energy_shift(data,x)*(1+F(x)**2)**2]
            for x in points]
    for k in range(5):
        pivot = next(i for i in range(k,5) if rows[i][k])
        rows[k],rows[pivot] = rows[pivot],rows[k]
        divisor = rows[k][k]
        rows[k] = [a/divisor for a in rows[k]]
        for i in range(5):
            if i!=k:
                scale = rows[i][k]
                rows[i] = [a-scale*b for a,b in zip(rows[i],rows[k])]
    return [rows[i][-1] for i in range(5)]


def family_bound(data,lower,gram_free):
    """Prove every real contact rotation's shift is at least lower."""
    lower,z = _exact(lower,10**9),_exact(gram_free,10**9)
    coefficients = shift_polynomial(data)
    p = [a-lower*b for a,b in zip(coefficients,(1,0,2,0,1))]
    matrix = [[p[0],p[1]/2,z],[p[1]/2,p[2]-2*z,p[3]/2],[z,p[3]/2,p[4]]]
    scale = lcm(*(a.denominator for row in matrix for a in row))
    psd = integer_psd([[int(a*scale) for a in row] for row in matrix])
    # [1,x,x²]^T Q [1,x,x²] exactly recovers the polynomial coefficients.
    recovered = [matrix[0][0],2*matrix[0][1],2*matrix[0][2]+matrix[1][1],
                 2*matrix[1][2],matrix[2][2]]
    if recovered!=p:raise ValueError('Exact quartic SOS identity failed')
    return {'accepted':True,'minimum_shift_lower':str(lower),'gram_free':str(z),
            'numerator_coefficients':list(map(str,coefficients)),'gram_psd':psd,
            'scope':'Restricted-family lower limit for the contact-rotation energy shift, including the x=infinity limit. It is not a lower bound on the unrestricted ground energy.'}


def replay(c,h,x,sites,remainder_cert=None):
    if type(sites) is not int or sites%2 or not 8<=sites<=10**9:
        raise ValueError('Even chain size between8 and10^9 required')
    x = _exact(x,10**6)
    q,r = divmod(sites,8)
    if not r and remainder_cert is not None:raise ValueError('No remainder is used')
    remainder = remainder_upper(remainder_cert,r) if r else None
    data = edge_reductions(c,h)
    shift = energy_shift(data,x)
    plain = q*F(data['energy8'])+(F(remainder['upper']) if remainder else 0)
    upper = plain+(q-1)*shift
    return {'accepted':True,'sites':sites,'blocks':q,'remainder':r,'x':str(x),
            'gate_energy_shift':str(shift),'upper':str(upper),'upper_per_site':str(upper/sites),
            'unfiltered_upper':str(plain),'improvement':str(plain-upper),'norm_factor':'1',
            'edge_reductions':data,'remainder_upper':remainder,
            'scope':'Physical depth-one contact-unitary circuit on fixed H8 building blocks. Each gate changes only a four-site neighborhood; other gates are disjoint from it, so energy increments add exactly. Endpoint two-site RDMs and the CAR contraction are freshly recomputed. No ground-state or circuit-family optimality claim.'}
