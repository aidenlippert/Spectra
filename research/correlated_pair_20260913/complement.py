"""Correlated-state residual and a restricted local-parent obstruction."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import decode,product
from research.correlated_pair_20260913.mps_exact import State,I
from research.molecular_collective_20260913.core import digest


def ldlt_positive(matrix):
    n=len(matrix);L=[[F(i==j) for j in range(n)] for i in range(n)];D=[]
    for i in range(n):
        d=matrix[i][i]-sum(L[i][k]**2*D[k] for k in range(i))
        if d<=0:return False,D+[d]
        D.append(d)
        for j in range(i+1,n):L[j][i]=(matrix[j][i]-sum(L[j][k]*L[i][k]*D[k] for k in range(i)))/d
    return True,D


def local_parent_obstruction(state):
    rows=[]
    for start in range(0,state.m-3,2):
        sites=list(range(start,start+4));groups={}
        for x in range(16):
            q=(sum((x>>i)&1 for i in (0,2)),sum((x>>i)&1 for i in (1,3)))
            groups.setdefault(q,[]).append(x)
        pivots=[];ok=True
        for basis in groups.values():
            matrix=[]
            for x in basis:
                row=[]
                for y in basis:
                    ops=[I]*state.m
                    for bit,site in enumerate(sites):
                        O=[0]*4;O[2*((y>>bit)&1)+((x>>bit)&1)]=1;ops[site]=tuple(O)
                    row.append(F(state.contract_factors(ops),state.norm_integer))
                matrix.append(row)
            if any(matrix[i][j]!=matrix[j][i] for i in range(len(basis)) for j in range(len(basis))):raise ValueError('RDM Hermiticity')
            good,ps=ldlt_positive(matrix);ok &= good;pivots.extend(ps)
        rows.append({'spatial_orbitals':[start//2,start//2+1],'RDM_dimension':16,'largest_charge_block':4,
                     'strictly_positive_definite':ok,'minimum_LDL_pivot':str(min(pivots))})
    return {'clusters':rows,'all_adjacent_two_orbital_RDMs_positive':all(x['strictly_positive_definite'] for x in rows),
            'conclusion':'If positive, no nonzero operator on any listed cluster exactly annihilates this rounded MPS. This obstructs exact frustration-free local parent factors in this support family only.',
            'complement_gamma_certified':False}


def run(data,cert,out,variance=True):
    start=time.monotonic();state=State(data,cert);h=decode(data['hamiltonian'],state.m,4);U=state.expectation(h)
    local=local_parent_obstruction(state);result={'fixture_sha256':digest(data),'state_sha256':digest(cert),'upper_Ha':str(U),'local_parent':local}
    out.write_text(json.dumps(result,indent=2)+'\n');print('Local parent check saved',flush=True)
    if variance:
        built=time.monotonic();h2=product(h,h);print(json.dumps({'H_squared_terms':len(h2),'H_squared_construction_seconds':time.monotonic()-built}),flush=True)
        variance_exact=state.expectation(h2)-U*U
        if variance_exact<0:raise AssertionError('Negative exact variance')
        result['residual']={'squared_norm_Ha2':str(variance_exact),'squared_norm_float_Ha2':float(variance_exact),'H_squared_terms':len(h2),'H_squared_max_degree':max(map(len,h2)),
            'claim':'Exact residual of normalized rounded MPS. No ground-state conclusion without a complement bound.'}
    result.update(wall_seconds=time.monotonic()-start,stats=state.stats)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'variance':result.get('residual',{}).get('squared_norm_float_Ha2'),'seconds':result['wall_seconds'],'local_obstruction':local['all_adjacent_two_orbital_RDMs_positive']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('output',type=Path);p.add_argument('--local-only',action='store_true');a=p.parse_args()
    run(json.loads(Path(a.fixture).read_text()),json.loads(Path(a.state).read_text()),a.output,not a.local_only)
