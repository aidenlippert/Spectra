"""Local Hubbard data and identity-seed MPS; no determinant enumeration."""
from fractions import Fraction as F
from math import comb
from research.direct_control_20260916.tensor_operator import build


def fixture(rungs=4):
    if type(rungs)!=int or not 2<=rungs<=8:raise ValueError('Bounded ladder size')
    sites=2*rungs;terms=[]
    def add(c,word):terms.append({'coefficient':str(c),'word':[list(x) for x in word]})
    edges=[(2*r,2*r+1) for r in range(rungs)]
    edges += [(2*r+s,2*(r+1)+s) for r in range(rungs-1) for s in (0,1)]
    for a,b in edges:
        for spin in (0,1):
            for i,j in ((a,b),(b,a)):add(-1,((1,2*i+spin),(0,2*j+spin)))
    for i in range(sites):add(-8,((1,2*i),(0,2*i),(1,2*i+1),(0,2*i+1)))
    add(8*rungs,())
    return {'modes':2*sites,'particles':2*rungs,'hamiltonian':terms,
            'description':'Particle-hole transformed repulsive half-filled bipartite Hubbard',
            'rungs':rungs,'spin_counts':[rungs,rungs],'t':'1','U':'8'}


def mpo(rungs=4):return build(fixture(rungs))


def numerical_mpo(op):
    import numpy as np
    arrays=[];charges=[[(0,0)]]
    for i,layer in enumerate(op['layers']):
        w=np.zeros((op['widths'][i],op['widths'][i+1],2,2));new=[None]*op['widths'][i+1]
        for l,r,mat,c in layer:
            for k,value in enumerate(mat):
                if not value:continue
                bra,ket=divmod(k,2);q=list(charges[-1][l]);q[i%2]+=bra-ket;q=tuple(q)
                if new[r] is not None and new[r]!=q:raise ValueError('Mixed MPO charge channel')
                new[r]=q;w[l,r,bra,ket]+=float(F(c)*value)
        if any(q is None for q in new):raise ValueError('Unreachable MPO bond')
        arrays.append(w);charges.append(new)
    if charges[-1]!=[(0,0)]:raise ValueError('Nonconserving MPO')
    return arrays,charges


def trace_seed(rungs=4):
    """Unnormalized identity in matrix space. Ground overlap is at least 1."""
    import numpy as np
    m=4*rungs;target=(rungs,rungs);reachable=[{(0,0)}];layers=[]
    for i in range(m):
        edges=[];following=set()
        for q in sorted(reachable[-1]):
            for s in (0,1):
                if i%2 and s!=q[0]-q[1]:continue
                qq=list(q);qq[i%2]+=s;qq=tuple(qq)
                if any(qq[j]>target[j] for j in (0,1)):continue
                edges.append((q,s,qq));following.add(qq)
        layers.append(edges);reachable.append(following)
    allowed={target};kept=[None]*m;charges=[None]*(m+1);charges[-1]=[target]
    for i in range(m-1,-1,-1):
        kept[i]=[e for e in layers[i] if e[2] in allowed]
        allowed={e[0] for e in kept[i]};charges[i]=sorted(allowed)
    arrays=[]
    for i,edges in enumerate(kept):
        left={q:j for j,q in enumerate(charges[i])};right={q:j for j,q in enumerate(charges[i+1])}
        a=np.zeros((len(left),2,len(right)))
        for q,s,r in edges:a[left[q],s,right[r]]=1.
        arrays.append(a)
    return arrays,charges


def seed_squared_norm(rungs=4):return comb(2*rungs,rungs)


def exact_seed_cert(rungs=4,denominator=1):
    """Standard-library charge-counter identity seed, for independent replay."""
    from research.molecular_collective_20260913.core import digest
    if type(denominator)!=int or denominator<1:raise ValueError('Positive denominator')
    data=fixture(rungs);m=data['modes'];target=(rungs,rungs);reachable={(0,0)};layers=[]
    for i in range(m):
        following=set();edges=[]
        for q in sorted(reachable):
            for s in (0,1):
                if i%2 and s!=q[0]-q[1]:continue
                qq=list(q);qq[i%2]+=s;qq=tuple(qq)
                if qq[0]>rungs or qq[1]>rungs:continue
                edges.append((q,s,qq));following.add(qq)
        layers.append(edges);reachable=following
    allowed={target};kept=[None]*m;charges=[None]*(m+1);charges[-1]=[target]
    for i in range(m-1,-1,-1):
        kept[i]=[e for e in layers[i] if e[2] in allowed]
        allowed={e[0] for e in kept[i]};charges[i]=sorted(allowed)
    tensors=[]
    for i,edges in enumerate(kept):
        left={q:j for j,q in enumerate(charges[i])};right={q:j for j,q in enumerate(charges[i+1])}
        tensors.append([[left[q],s,right[r],denominator] for q,s,r in edges])
    return {'kind':'integer_charge_mps_v1','fixture_sha256':digest(data),'modes':m,
      'particles':data['particles'],'spin_counts':data['spin_counts'],'denominator':denominator,
      'bond_charges':[[list(q) for q in layer] for layer in charges],'tensors':tensors}
