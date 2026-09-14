"""Exact extensive upper states from number-conserving boundary filters.

Each filter joins a previously correlated cluster to a fresh H8 block. Endpoint
observables persist under the merge, allowing an exact additive energy formula.
Only the original bounded H8 orbit state is reconstructed; no global state or
sixteen-site sector is expanded.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm

from experiments.marginal_local_hubbard_block import _exact, _terms
from experiments.marginal_symbolic import canonical, mono, product, add, scale
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
from experiments.marginal_tiled_upper import _uniform_h, replay_upper8, remainder_upper
from experiments.marginal_transfer_verify import apply_word
from results.marginal_graded_hubbard8.discovery.singlet_moment_energy import upper_monomials


def verify_edge_identity(sites,edge,spin):
    """Exactly check c H c†+c† H c=H+T_sigma+4(n_opposite-2D)."""
    if (type(sites) is not int or sites not in (2,4,8)
            or type(edge) is not int or edge not in (0,sites-1)
            or type(spin) is not int or spin not in (0,1)):
        raise ValueError('Bounded cluster endpoint and spin required')
    mode,other = 2*edge+spin,2*edge+1-spin
    neighbor = 2*(1 if edge==0 else sites-2)+spin
    c,cd = mono(((0,mode),)),mono(((1,mode),))
    h = canonical(dict(_terms(sites,F(4),F(1))))
    actual = add(product(product(c,h),cd),product(product(cd,h),c))
    n = mono(((1,mode),(0,mode)))
    opposite = mono(((1,other),(0,other)))
    hopping = add(mono(((1,mode),(0,neighbor))),mono(((1,neighbor),(0,mode))))
    expected = add(h,hopping,scale(opposite,4),scale(product(n,opposite),-8))
    residual = canonical(add(actual,scale(expected,-1)))
    if residual:
        raise ValueError('Charged endpoint CAR identity failed')
    return {'sites':sites,'edge':edge,'spin':spin,'residual_terms':0}


def merge_shift(delta,eta):
    # delta is a derived exact observable, not a user-controlled certificate.
    if type(delta) not in (int,str,F):
        raise ValueError('Exact boundary energy increment required')
    delta = F(delta)
    eta = _exact(eta,10**6)
    return (-2*eta+delta*eta*eta)/(1+eta*eta)


def boundary_data(c,h):
    checked = replay_upper8(c,h)
    oracle = SymmetryMomentOracle(_uniform_h(8))
    vectors = [oracle.compress({int(s):a for s,a in v.items()})
               for v in c['embedding']['basis']]
    table = upper_monomials(c['upper_chebyshev_coefficients'])
    state = {}
    for row in reversed(table):
        state = oracle.apply(state)
        for a,v in zip(row,vectors):
            for s,b in v.items():
                state[s] = state.get(s,0)+a*b
        state = {s:a for s,a in state.items() if a}
        if len(state) > 4096:
            raise ValueError('Derived orbit support exceeds existing sparse budget')
    if not state:
        raise ValueError('Zero polynomial trial state')
    denominator = lcm(*(a.denominator for a in state.values()))
    integers = {s:int(a*denominator) for s,a in state.items()}
    divisor = reduce(gcd,(abs(a) for a in integers.values()))
    state = {s:a//divisor for s,a in integers.items()}
    norm = oracle.dot(state,state)
    energy = F(oracle.dot(state,oracle.apply(state)),norm)
    if energy != F(checked['upper']):
        raise ValueError('Horner physical state disagrees with independent moment quotient')

    def amplitude(s):
        entry = oracle.orbit(s)
        return 0 if entry is None else state.get(entry[0],0)*entry[1]

    doublons,hopping,occupations = [0,0],[0,0],[[0,0],[0,0]]
    streamed = 0
    for r,a in state.items():
        # Stream the at-most-four members of each orbit. Never construct a
        # full determinant vector or a charged-sector wavefunction.
        for s,phase in oracle.orbit(r)[3].items():
            streamed += 1
            amp = a*phase
            for j,i in enumerate((0,7)):
                doublons[j] += amp*amp*int(((s>>(2*i))&3)==3)
                for spin in (0,1):
                    occupations[j][spin] += amp*amp*((s>>(2*i+spin))&1)
            for j,i in enumerate((0,6)):
                for spin in (0,1):
                    x,y = 2*i+spin,2*(i+1)+spin
                    for word in (((1,x),(0,y)),((1,y),(0,x))):
                        result = apply_word(word,s)
                        if result:
                            hopping[j] += amp*result[1]*amplitude(result[0])
    doublons = [F(x,norm) for x in doublons]
    hopping = [F(x,norm) for x in hopping]
    occupations = [[F(x,norm) for x in row] for row in occupations]
    if occupations != [[F(1,2)]*2]*2:
        raise ValueError('Each endpoint spin must have exact occupation1/2')
    identities = [verify_edge_identity(8,edge,spin) for edge in (0,7) for spin in (0,1)]
    delta = sum(hopping)/2+4*(1-2*sum(doublons))
    return {'accepted':True,'energy8':str(energy),'norm':str(norm),
            'edge_doublons':list(map(str,doublons)),
            'edge_hopping':list(map(str,hopping)),
            'edge_spin_occupations':[[str(x) for x in row] for row in occupations],
            'delta':str(delta),'car_identities':identities,
            'orbit_support':len(state),'streamed_determinants':streamed,
            'unique_source_actions':len(oracle.base.cache),'upper_replay':checked,
            'scope':'Exact endpoint observables of the validated physical H8 polynomial state, checked against the independently recomputed moment energy. No reflection or eigenstate assumption.'}


def replay(c,h,eta,sites,remainder_cert=None):
    if type(sites) is not int or sites%2 or not 8 <= sites <= 10**9:
        raise ValueError('Even chain length between8 and10^9 required')
    eta = _exact(eta,10**6)
    q,r = divmod(sites,8)
    if not r and remainder_cert is not None:
        raise ValueError('No remainder is used for a multiple of8')
    remainder = remainder_upper(remainder_cert,r) if r else None
    data = boundary_data(c,h)
    shift = merge_shift(data['delta'],eta)
    unfiltered = q*F(data['energy8'])+(F(remainder['upper']) if remainder else 0)
    upper = unfiltered+(q-1)*shift
    return {'accepted':True,'sites':sites,'blocks':q,'remainder':r,
            'eta':str(eta),'merge_energy_shift':str(shift),
            'upper':str(upper),'upper_per_site':str(upper/sites),
            'unfiltered_upper':str(unfiltered),'improvement':str(unfiltered-upper),
            'filter_norm_base':str(1+eta*eta),'filter_norm_exponent':q-1,
            'boundary_data':data,'remainder_upper':remainder,
            'scope':'Physical correlated upper state: apply I-eta*h at every cut between consecutive H8 blocks, then append an unfiltered fixed-particle remainder if present. The merged cluster preserves total spin populations; remote endpoint observables persist exactly. Each merge adds the verified constant shift. Norm and state remain compact recipes; no global vector or exponentially large norm integer is expanded. This is not a claim of efficient physical postselection or state preparation.'}
