"""Exact PH density allocation obstruction and paired density-conditioned CAR.

This supersedes the preliminary maxflow.py/closure.py experiments. A cut
certifies failure of the attractive-density parent allocation family only.
"""
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from experiments.marginal_symbolic import decode,mono,canonical,product,add,scale,adj
from research.correlated_pair_20260913.self_consistent.fixed_guide import guide,ph
from research.molecular_collective_20260913.core import digest


def number(i):return mono(((1,i),(0,i)))


def density_inputs(data):
    mu,weights,V,diag=guide(data);hp=ph(decode(data['hamiltonian'],data['modes'],4),data['particles']);pairs={}
    for i,j in combinations(range(data['modes']),2):
        word,sign=next(iter(product(number(i),number(j)).items()))
        value=hp.get(word,F(0))/sign
        if value:pairs[i,j]=value
    return mu,weights,pairs,diag


def allocation(weights,pairs):
    if any(w<0 for w in weights):raise ValueError('Nonnegative capacities required')
    edges=[(i,j,-u) for (i,j),u in sorted(pairs.items()) if u<0];need=sum((c for i,j,c in edges),F(0));big=need+1
    k=len(edges);m=len(weights);sink=k+m+1;cap={};adjacency={i:set() for i in range(sink+1)}
    def edge(u,v,c):cap[u,v]=c;cap[v,u]=F(0);adjacency[u].add(v);adjacency[v].add(u)
    for a,(i,j,c) in enumerate(edges,1):edge(0,a,c);edge(a,k+1+i,big);edge(a,k+1+j,big)
    for i,w in enumerate(weights):edge(k+1+i,sink,w)
    rem=cap.copy();got=F(0)
    while True:
        parent={0:None};queue=[0]
        for u in queue:
            for v in sorted(adjacency[u]):
                if v not in parent and rem[u,v]>0:parent[v]=u;queue.append(v)
        if sink not in parent:break
        delta=need+1;v=sink
        while v:delta=min(delta,rem[parent[v],v]);v=parent[v]
        v=sink
        while v:u=parent[v];rem[u,v]-=delta;rem[v,u]+=delta;v=u
        got+=delta
    reach=set(parent);cut=sum((c for (u,v),c in cap.items() if u in reach and v not in reach),F(0))
    if cut!=got:raise AssertionError('Flow and cut disagree')
    # Extract a simple independently verifiable subset inequality. For any S,
    # sum_{i,j in S}(-Uij) <= sum_{i in S}ei is necessary for allocations.
    S=sorted(i for i in range(m) if k+1+i in reach)
    required=sum((c for i,j,c in edges if i in S and j in S),F(0));capacity=sum((weights[i] for i in S),F(0))
    if got<need and required<=capacity:raise AssertionError('Missing exact subset obstruction')
    alloc=[{'pair':[i,j],'to_first':str(big-rem[a,k+1+i]),'to_second':str(big-rem[a,k+1+j]),'required':str(c)} for a,(i,j,c) in enumerate(edges,1)]
    if got==need:
        for row in alloc:
            if F(row['to_first'])+F(row['to_second'])!=F(row['required']):raise AssertionError('Unfilled pair')
        for i,w in enumerate(weights):
            used=sum((F(row['to_first'] if row['pair'][0]==i else row['to_second']) for row in alloc if i in row['pair']),F(0))
            if used>w:raise AssertionError('Overfilled node')
    return {'feasible':got==need,'flow_Ha':str(got),'required_Ha':str(need),'cut_capacity_Ha':str(cut),
            'subset':S,'subset_required_Ha':str(required),'subset_capacity_Ha':str(capacity),
            'subset_deficit_Ha':str(required-capacity),'allocation':alloc if got==need else None}


def paired(B,theta):
    f=add(B,theta)
    return canonical(add(product(adj(f),f),product(theta,adj(theta))))


def closure_demo():
    # Allowed quintic correction: n5 times a cubic operator commuting with n5.
    nj=number(5);nu=mono(((1,0),(0,2),(0,1)),F(1,7))
    if product(nj,nu)!=product(nu,nj):raise AssertionError('Conditioner commutator')
    B=add(mono(((0,0),)),scale(product(nj,mono(((0,1),))),F(1,3)))
    theta=add(mono(((1,1),(1,2),(1,3)),F(1,5)),mono(((0,4),),F(1,11)),product(nj,nu))
    generated=paired(B,theta);six={w:c for w,c in generated.items() if len(w)==6}
    if not six or max(map(len,generated))!=6:raise AssertionError('Expected retained sixth degree closure')
    # The four-mode strong-density example: U=2d1+d1^2+2d2+d2^2.
    eps=F(1,4);d=F(1,2);U=4*d+2*d*d;den=4+U
    base=add(*(number(i) for i in range(4)),scale(product(number(0),number(1)),U))
    v=mono(tuple((1,i) for i in range(4)),eps);target=add(base,v,adj(v));sos={}
    for i in range(4):
        B=mono(((0,i),));mult=F(1)
        if i<2:B=add(B,scale(product(number(1-i),B),d));mult+=d
        theta=mono(tuple((1,j) for j in range(4) if j!=i),((-1)**i)*eps*mult/den)
        sos=add(sos,paired(B,theta))
    rem=canonical(add(target,scale(sos,-1)));b=rem.pop((),F(0));err=sum(map(abs,rem.values()),F(0))
    return {'six_mode_allowed_closure':{'terms':len(generated),'degree6_terms':len(six),'max_degree':6,
                'quartic_truncation_L1_operator_penalty':str(sum(map(abs,six.values()),F(0)))},
            'four_mode_resonant_example':{'U':str(U),'epsilon':str(eps),'positive_factors':8,'lower_Ha':str(b-err),
                'base_Ha':str(b),'residual_L1_Ha':str(err),'full_CAR_residual_terms':len(rem),
                'exact_ground_energy_formula':'((4+U)-sqrt((4+U)^2+4*epsilon^2))/2; all other diagonal sectors nonnegative'},
            'scope':'Exact toy mechanism and closure, not a molecular ground-energy certificate.'}


def run():
    out=Path('results/correlated_pair_20260913/density_conditioned');out.mkdir(parents=True,exist_ok=True);report={'closure':closure_demo(),'cases':{}}
    for case in ('h6','h8'):
        data=json.loads(Path(f'results/certificate_scaling/active_space_ladder/{case}/fixture.json').read_text());mu,e,U,diag=density_inputs(data)
        report['cases'][case]={'fixture_sha256':digest(data),'chemical_potential':str(mu),'weights':list(map(str,e)),
                'pair_density_terms':len(U),'attractive_pairs':sum(u<0 for u in U.values()),'flow':allocation(e,U),
                'scope':'Restricted PH attractive-density parent family; positive U dropped, remaining offdiagonal H not tested by this cut.'}
    (out/'exact_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':run()
