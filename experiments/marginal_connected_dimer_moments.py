"""Exact fixed-order moments from connected intervals of Hubbard dimers.

The normalized boundary is a product of adjacent two-site valence singlets.
Formal log moments are additive when an inter-dimer bond is removed. A
connected interval of b dimers therefore uses all b-1 cuts. Returning to the
boundary charge requires an even number of hops across each cut, so its
order-n cumulant weight vanishes for n<2(b-1). At order8, clusters of at most
10 sites suffice for any even chain length. This does not determine the
spectrum or ground energy from finitely many moments.
"""
from fractions import Fraction as F
from math import comb
from experiments.marginal_clifford_moments import CliffordMomentOracle
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_symmetry_moments import simple


def cumulants(moments):
    if not moments or moments[0]!=1:raise ValueError('Normalized moment zero must equal one')
    result=[0]
    for n in range(1,len(moments)):
        result.append(simple(moments[n]-sum(comb(n-1,j-1)*result[j]*moments[n-j] for j in range(1,n))))
    return result


def from_cumulants(values):
    result=[1]
    for n in range(1,len(values)):
        result.append(simple(sum(comb(n-1,j-1)*values[j]*result[n-j] for j in range(1,n+1))))
    return result


def compile_moments(sites,U,t,maximum_order=8):
    if type(sites) is not int or sites<2 or sites>10**9 or sites%2:
        raise ValueError('Even site count in2..10^9 required')
    if type(maximum_order) is not int or not 0<=maximum_order<=8:
        raise ValueError('Connected moment order must lie in0..8')
    for value in (U,t):
        if type(value) not in (str,int,F):raise ValueError('Exact Hubbard parameters required')
        q=F(value)
        if abs(q)>10**6 or q.denominator>10**6:raise ValueError('Hubbard parameter bound exceeded')
    U,t=F(U),F(t)
    if U<0 or t<0:raise ValueError('Nonnegative U and t required for this declared model family')
    dimers=sites//2;limit=min(dimers,maximum_order//2+1) if maximum_order else 0
    kappas=[[0]*(maximum_order+1)];weights=[];work=[]
    for b in range(1,limit+1):
        oracle=CliffordMomentOracle(build(2*b,U,t));matrices,receipt=oracle.normalized_dimer_moments(maximum_order)
        kappa=cumulants([m[0][0] for m in matrices]);kappas.append(kappa)
        weight=[kappa[n]-2*kappas[b-1][n]+(kappas[b-2][n] if b>=2 else 0) for n in range(maximum_order+1)]
        if any(weight[n] for n in range(min(maximum_order+1,2*(b-1)))):
            raise ValueError('Connected interval failed the charge-return order gate')
        weights.append(weight);work.append(dict(receipt,cluster_dimers=b))
    total=[0]*(maximum_order+1)
    for n in range(1,maximum_order+1):
        total[n]=simple(sum((dimers-b+1)*weights[b-1][n] for b in range(1,min(limit,n//2+1)+1)))
    moments=from_cumulants(total)
    return moments,{'sites':sites,'dimers':dimers,'U':str(U),'t':str(t),'maximum_order':maximum_order,
        'normalized_boundary':True,'maximum_cluster_sites':2*limit,'cluster_count':limit,
        'total_source_frame_states':sum(w['source_frame_states'] for w in work),
        'largest_cluster_source_states':max((w['source_frame_states'] for w in work),default=0),
        'cluster_work':work,'connected_cumulant_weights':[[str(x) for x in row] for row in weights],
        'chain_cumulants':[str(x) for x in total],
        'scope':'Exact normalized moments of the uniform open half-filled Hubbard chain in one product-of-valence-singlets boundary, through the declared fixed order. All interval weights are recomputed on at most10sites. No full chain configurations or Hamiltonian matrix; no ground-energy, general-boundary or universal representability claim.'}
