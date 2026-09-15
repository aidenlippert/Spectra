"""Physical block-state adaptation for the complete filtered-chain objective."""
from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce

from experiments.marginal_boundary_unitary import _physical_source, _state_reductions, model_shift
from experiments.marginal_local_hubbard_block import _exact
from experiments.marginal_symbolic import decode,encode,mono,product,add,scale
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
from experiments.marginal_tiled_upper import _uniform_h
from experiments.marginal_polynomial_sos import integer_psd


def effective_model(eta):
    eta=_exact(eta)
    w=eta*eta/(1+eta*eta)
    h=decode(_uniform_h(8)['hamiltonian'],16,4)
    boundary={}
    for edge,neighbor in ((0,1),(7,6)):
        up=mono(((1,2*edge),(0,2*edge)))
        down=mono(((1,2*edge+1),(0,2*edge+1)))
        boundary=add(boundary,scale(product(up,down),-8),scale(up,4),scale(down,4),mono((),-4))
        for spin in (0,1):
            a,b=2*edge+spin,2*neighbor+spin
            boundary=add(boundary,mono(((1,a),(0,b)),F(1,2)),mono(((1,b),(0,a)),F(1,2)))
    model={'modes':16,'particles':8,'hamiltonian':encode(add(h,scale(boundary,w)))}
    return model,w


def _primitive(state):
    state={s:a for s,a in state.items() if a}
    if not state or len(state)>4096:raise ValueError('Nonzero bounded physical orbit state required')
    denominator=lcm(*(a.denominator for a in state.values()))
    integers={s:int(a*denominator) for s,a in state.items()}
    divisor=reduce(gcd,(abs(a) for a in integers.values()))
    return {s:a//divisor for s,a in integers.items()}


def krylov_lower(c,h,eta,degree,lower):
    """Prove a lower Rayleigh limit in the specified finite physical span only."""
    if type(degree) is not int or not 0<=degree<=8:
        raise ValueError('Krylov degree must be between0 and8')
    eta=_exact(eta);lower=_exact(lower,10**12)
    _,source,source_receipt=_physical_source(c,h)
    model,w=effective_model(eta)
    oracle=SymmetryMomentOracle(model)
    vectors=[source]
    images=[]
    for j in range(degree+1):
        image=oracle.apply(vectors[-1]);images.append(image)
        if j<degree:
            shifted={s:image.get(s,0)+4*vectors[-1].get(s,0)
                     for s in image.keys()|vectors[-1].keys()}
            vectors.append(_primitive(shifted))
    # Each vector is rescaled independently by a nonzero rational. This
    # preserves the polynomial span and controls intermediate denominators.
    gram=[[oracle.dot(a,b) for b in vectors] for a in vectors]
    energy=[[oracle.dot(a,b) for b in images] for a in vectors]
    size=len(vectors)
    if any(energy[i][j]!=energy[j][i] for i in range(size) for j in range(size)):
        raise ValueError('Krylov energy Hermiticity failed')
    shifted=[[energy[i][j]-lower*gram[i][j] for j in range(size)] for i in range(size)]
    denominator=lcm(*(a.denominator for row in shifted for a in row))
    psd=integer_psd([[int(a*denominator) for a in row] for row in shifted])
    return {'accepted':True,'eta':str(eta),'degree':degree,'effective_energy_lower':str(lower),
            'thermodynamic_density_lower':str((lower+4*w-2*eta/(1+eta*eta))/8),
            'gram_dimension':size,'matrix_psd':psd,'source_upper_replay':source_receipt,
            'scope':'Lower bound only on K_eta Rayleigh quotients in span{(K_eta+4I)^j source:0<=j<=degree}, at the specified eta. Not a ground-energy lower bound, not a limit across eta or arbitrary block states.'}


def replay(c,h,eta,coefficients,sites):
    """Apply a real polynomial in K_eta+4I to the validated physical source."""
    if type(sites) is not int or sites%8 or not 8<=sites<=10**9:
        raise ValueError('Adapted chain currently requires multiples of8 sites')
    eta=_exact(eta)
    if type(coefficients) is not list or not 1<=len(coefficients)<=9:
        raise ValueError('Polynomial degree must be between0 and8')
    coefficients=[_exact(a,10**12) for a in coefficients]
    if not any(coefficients):raise ValueError('Nonzero polynomial recipe required')
    original,source,source_receipt=_physical_source(c,h)
    model,w=effective_model(eta)
    effective=SymmetryMomentOracle(model)
    # Orbit conventions are identical; the new oracle independently verifies
    # the actual effective Hamiltonian's conservation and discrete symmetries.
    state={}
    maximum_support=0
    for coefficient in reversed(coefficients):
        moved=effective.apply(state)
        state={s:moved.get(s,0)+4*state.get(s,0)+coefficient*source.get(s,0)
               for s in moved.keys()|state.keys()|source.keys()}
        state={s:a for s,a in state.items() if a}
        maximum_support=max(maximum_support,len(state))
        if maximum_support>4096:raise ValueError('Existing orbit support budget exceeded')
    state=_primitive(state)
    data=_state_reductions(original,state)
    norm=original.dot(state,state)
    k_energy=F(effective.dot(state,effective.apply(state)),norm)
    for name,site in (('right',1),('left',0)):
        for spin in (0,1):
            value=sum(F(data[name][i][i])*((i>>(2*site+spin))&1) for i in range(16))
            if value!=F(1,2):raise ValueError('Adapted endpoint half occupation failed')
    e=F(data['energy8'])
    # This identity yields the thermodynamic objective, not finite-size
    # optimality. Finite-chain assembly still keeps exactly q-1 cuts.
    objective=k_energy+4*w-2*eta/(1+eta*eta)
    shift=objective-e
    direct=model_shift(data,eta,'linear_filter')
    if shift!=direct:raise ValueError('Effective objective disagrees with fresh four-site CAR contraction')
    q=sites//8
    upper=q*e+(q-1)*shift
    return {'accepted':True,'sites':sites,'blocks':q,'eta':str(eta),
            'polynomial_coefficients':list(map(str,coefficients)),
            'polynomial_variable':'K_eta+4I','degree':len(coefficients)-1,
            'block_energy':str(e),'effective_energy':str(k_energy),
            'thermodynamic_objective_per_block':str(objective),
            'merge_shift':str(shift),'upper':str(upper),'upper_per_site':str(upper/sites),
            'edge_reductions':data,'source_upper_replay':source_receipt,
            'effective_symmetry':effective.symmetry_receipt,
            'maximum_orbit_support':maximum_support,'source_actions':len(effective.base.cache),
            'scope':'Exact physical adapted block from a bounded degree<=8 polynomial in the verified effective Hamiltonian. Fresh original-H energy, endpoint RDMs, and independent four-site filter contraction agree. This is a variational upper certificate, not proof of global or block-family optimality.'}
