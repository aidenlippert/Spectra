"""Exact rational feasibility for pair-density parent factors.

For every attractive density edge U_ij<0 allocate c_ij=-U_ij between its two
endpoints.  The allocation is feasible when node capacities e_i can absorb
all allocations.  This is a finite max-flow certificate, independent of any
many-body basis.  Positive U terms are separate positive squares.
"""
from fractions import Fraction as F
raise RuntimeError('Rejected exploratory prototype: use density_conditioned.exact_density for correct PH signs, capacities and cut verification.')
from math import gcd

def _lcm(a,b): return a//gcd(a,b)*b

def _flow(cap, source, sink):
    n=len(cap); f=[[0]*n for _ in range(n)]
    total=0
    while True:
        parent=[-1]*n; parent[source]=source; q=[source]
        for u in q:
            for v,c in enumerate(cap[u]):
                if parent[v]<0 and c-f[u][v]>0: parent[v]=u;q.append(v)
        if parent[sink]<0: break
        d=10**100;v=sink
        while v!=source:u=parent[v];d=min(d,cap[u][v]-f[u][v]);v=u
        v=sink
        while v!=source:u=parent[v];f[u][v]+=d;f[v][u]-=d;v=u
        total+=d
    return total,f

def feasible(e, negative_edges):
    e=list(map(F,e)); edges=[(i,j,F(c)) for i,j,c in negative_edges if F(c)>0]
    den=1
    for x in e+[c for _,_,c in edges]: den=_lcm(den,x.denominator)
    E=[int(x*den) for x in e]; C=[int(c*den) for _,_,c in edges]
    m=len(e); source=0; edge0=1; node0=edge0+len(edges); sink=node0+m; n=sink+1
    cap=[[0]*n for _ in range(n)]
    for k,c in enumerate(C):
        cap[source][edge0+k]=c; i,j,_=edges[k]; cap[edge0+k][node0+i]=10**30;cap[edge0+k][node0+j]=10**30
    if any(c<0 for c in E): raise ValueError('negative node capacity')
    for i,c in enumerate(E): cap[node0+i][sink]=c
    got,flow_matrix=_flow(cap,source,sink); need=sum(C)
    reach={source}; changed=True
    while changed:
      changed=False
      for u in list(reach):
        for v,c in enumerate(cap[u]):
          if v not in reach and c-flow_matrix[u][v]>0: reach.add(v);changed=True
    return {'feasible':got==need,'flow':F(got,den),'required':F(need,den),'scale':den,
            'reachable_nodes':sorted(reach),'cut_capacity':F(sum(cap[u][v] for u in reach for v in range(n) if v not in reach),den),
            'edges':len(edges),'nodes':m}

def fixture_terms(data):
    one={}; pair={}
    for x in data['hamiltonian']:
        w=tuple(tuple(z) for z in x['word']); c=F(x['coefficient'])
        if len(w)==2 and w[0][0]==1 and w[1][0]==0 and w[0][1]==w[1][1]: one[w[0][1]]=one.get(w[0][1],F(0))+c
        if len(w)==4 and sorted(w)==sorted([(1,w[0][1]),(0,w[0][1]),(1,w[-1][1]),(0,w[-1][1])]) and w[0][1]!=w[-1][1]:
            i,j=sorted((w[0][1],w[-1][1]));pair[i,j]=pair.get((i,j),F(0))+c
    return one,pair

def particle_hole_density(data):
    """Transform the diagonal density polynomial around occupied 0..N-1."""
    one,pair=fixture_terms(data); n=data['particles']; m=data['modes']; occ=set(range(n))
    # n_i = 1-h_i on the reference-occupied modes, h_i otherwise.
    constant=sum(c for i,c in one.items() if i in occ)
    linear={}; U={}
    for i,c in one.items(): linear[i]=linear.get(i,F(0))+(-c if i in occ else c)
    for (i,j),c in pair.items():
        si=-1 if i in occ else 1; sj=-1 if j in occ else 1
        constant += c if i in occ and j in occ else F(0)
        linear[i]=linear.get(i,F(0))+(-c if i in occ else (c if j in occ else F(0)))
        linear[j]=linear.get(j,F(0))+(-c if j in occ else (c if i in occ else F(0)))
        U[i,j]=U.get((i,j),F(0))+si*sj*c
    return constant,linear,U

def shifted_particle_hole(data):
    c,e,U=particle_hole_density(data); mu=-min(e.values(),default=F(0))+F(1,10**8)
    occ=set(range(data['particles'])); lo=-min((v for i,v in e.items() if i not in occ),default=F(0)); hi=min((v for i,v in e.items() if i in occ),default=F(0))
    if lo>=hi: raise ValueError('no chemical-potential interval for positive particle-hole capacities')
    mu=(lo+hi)/2
    return c,{i:v-mu if i in occ else v+mu for i,v in e.items()},U,mu
