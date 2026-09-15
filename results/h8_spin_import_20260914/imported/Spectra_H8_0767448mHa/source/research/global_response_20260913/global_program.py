"""One simultaneous response on the union of two occupation sectors.

The accepting path uses orbital matrices and local CAR blocks, never a
particle-sector basis. Terminal positivity is a separate proof obligation.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys
import time

from experiments.marginal_symbolic import mono,product,add,scale
from research.compact_response_20260913 import program
from research.composable_response_20260913 import joint,spatial_gap,coupling
from research.molecular_collective_20260913.core import extract,tail_replay,digest
from research.positive_response_20260913.coercivity import spatial_inputs,spatial_norm
from research.positive_response_20260913.molecular_diagnostic import spatial_one_body

ROOT=program.ROOT
OUT=ROOT/'results/global_response_20260913'


def q_orbital(j):
    return product(mono(((1,2*j),(0,2*j))),mono(((1,2*j+1),(0,2*j+1))))


def second_bare_upper_and_coupling(data,tail,b):
    s=data['modes']//2;j=s-2;n=data['particles']-2;rem=[i for i in range(s) if i!=j]
    p=extract(data,tail['center_number']);c,t=spatial_one_body(p);Ls,weights,_=spatial_inputs(data,tail)
    upper=c+2*t[j][j]+program.one_body_endpoint([[t[i][k] for k in rem] for i in rem],n,True)
    for w,L in zip(weights,Ls):
        A=[[L[i][k] for k in rem] for i in rem];vv=sum(L[i][j]**2 for i in rem)
        upper+=w*((n*spatial_norm(A)+2*abs(L[j][j]))**2/2+vv)
    upper+=F(tail_replay(data,tail)['upper_operator_shift_Ha'])-b
    q=q_orbital(j);cross=product(product(q,p['h']),add(mono(()),scale(q,-1)))
    return upper,sum(abs(c) for c in cross.values()),len(cross)


def bounds(data,tail,b,first_gap,joint_gap):
    first=program.bounds(data,tail,b)
    if first_gap['orbital']!=data['modes']//2-1 or joint_gap['last_spatial_orbitals']!=2:
        raise ValueError('Wrong fixed occupation partition')
    delta1=F(spatial_gap.check(data,tail,first_gap)['lower_H_QQ_Ha'])-b
    delta12=F(joint.check_gap(data,tail,joint_gap)['lower_H_joint_Ha'])-b
    if delta12>0:
        # Q=Q1+R2 with R2=(I-Q1)Q2. Internal off-diagonal coupling
        # must be paid in the spectral upper endpoint, not dropped.
        M2,g2,terms=second_bare_upper_and_coupling(data,tail,b)
        cross=coupling.prove(data);nu=F(cross['coupling_norm_squared_Ha2'])
        M=(first['M']+M2+coupling.sqrt_up((first['M']-M2)**2+4*nu))/2
        # Q1HP and R2HP have orthogonal output ranges. P is a subspace
        # of both I-Q1 and I-Q2, so these external bounds restrict safely.
        g=coupling.sqrt_up(first['coupling_norm']**2+g2**2)
        return {'partition':'union_last_two_double','delta':delta12,'M':M,'g':g,
            'first_diagonal_M':first['M'],'second_diagonal_M':M2,'internal_nu':nu,
            'external_g1':first['coupling_norm'],'external_g2':g2,
            'bare_second_coupling_terms':terms,'local_coupling_work':cross}
    if delta1>0:
        return {'partition':'last_double','delta':delta1,'M':first['M'],'g':first['coupling_norm'],
            'joint_gap_inconclusive':delta12,'local_coupling_work':None}
    return {'partition':'none','delta':delta1,'M':first['M'],'g':first['coupling_norm'],
        'joint_gap_inconclusive':delta12,'local_coupling_work':None}


def scalar_program(bound,budget):
    if not F(0)<budget<=F(1,1000):raise ValueError('Positive bounded response error required')
    d=program.floor_grid(bound['delta']);M=program.ceil_grid(bound['M']);g=program.ceil_grid(bound['g'])
    if not 0<d<M or g<=0:return {'status':'no_certified_denominator','delta_Ha':str(d)}
    z=(M+d)/(M-d)
    for k in range(1,257):
        eta=g*g/(d*program.chebyshev(k,z)**2)
        if eta<=budget:
            return {'status':'certified_response','delta_Ha':str(d),'M_Ha':str(M),
                'coupling_norm_Ha':str(g),'order':k,'residual_penalty_Ha':str(eta),
                'residual_budget_Ha':str(budget)}
    return {'status':'order_cap_exhausted','order_cap':256,'delta_Ha':str(d),
        'M_Ha':str(M),'coupling_norm_Ha':str(g),'penalty_at_cap_Ha':str(eta)}


def check(data,tail,cert):
    start=time.monotonic()
    keys={'kind','fixture_sha256','tail_sha256','target_Ha','budget_Ha','first_sector','joint_sector','partition','response'}
    if set(cert)!=keys or cert['kind']!='global_occupation_response_v1':raise ValueError('Unknown global response schema')
    if cert['fixture_sha256']!=digest(data) or cert['tail_sha256']!=digest(tail):raise ValueError('Global input binding failed')
    if any(type(cert[k]) is not str for k in ('target_Ha','budget_Ha')):raise ValueError('Exact scalar inputs required')
    b=F(cert['target_Ha']);budget=F(cert['budget_Ha']);derived=bounds(data,tail,b,cert['first_sector'],cert['joint_sector'])
    expected=scalar_program(derived,budget)
    if cert['partition']!=derived['partition'] or cert['response']!=expected:raise ValueError('Global response did not reproduce')
    k=expected.get('order');s=data['modes']//2
    return {'status':expected['status'],'partition':cert['partition'],'target_Ha':str(b),
        'gap_lower_Ha':str(derived['delta']),'gap_lower_float_Ha':float(derived['delta']),
        'order':k,'H_actions_per_K_vector':2*k+1 if k else None,
        'D_actions_per_supplied_response_rhs':k-1 if k else None,
        'residual_penalty_Ha':expected.get('residual_penalty_Ha'),
        'many_body_states_enumerated':0,'many_body_matrix_entries':0,
        'new_orbital_matrix_dimension':s,'shared_tail_coefficient_dimension':s*(s+1)//2,
        'internal_coupling_squared_bound_Ha2':str(derived['internal_nu']) if 'internal_nu' in derived else None,
        'external_coupling_bounds_Ha':[str(derived[k]) for k in ('external_g1','external_g2')] if 'external_g1' in derived else None,
        'local_coupling_work':derived['local_coupling_work'],
        'terminal_positivity_proved_here':False,'replay_seconds':time.monotonic()-start}


def propose(data,tail,b,budget=F(1,10**6)):
    start=time.monotonic();fg,fs=spatial_gap.propose(data,tail,data['modes']//2-1);jg,js=joint.propose_gap(data,tail)
    derived=bounds(data,tail,b,fg,jg)
    cert={'kind':'global_occupation_response_v1','fixture_sha256':digest(data),'tail_sha256':digest(tail),
        'target_Ha':str(b),'budget_Ha':str(budget),'first_sector':fg,'joint_sector':jg,
        'partition':derived['partition'],'response':scalar_program(derived,budget)}
    return cert,{'first_gap_discovery':fs,'joint_gap_discovery':js,'construction_seconds':time.monotonic()-start}


def load_case(name):
    if name in ('h6','fresh_h6_1p6'):
        data,tail,_,first=joint.load_case(name);b=F(first['target_Ha'])
        return data,tail,b
    if name in ('h8','h8_legacy_target'):
        data,tail,_,first=joint.load_case('h8')
        if name=='h8_legacy_target':b=F(first['target_Ha'])
        else:
            ref=json.loads((ROOT/'results/certificate_scaling/cubic_precision/intervals/final_h8.json').read_text())
            b=F(ref['upper'])-F(1,1000)
        return data,tail,b
    if name=='fresh_h6_1p73':
        folder=OUT/name;data=json.loads((folder/'fixture.json').read_text());tail=json.loads((folder/'tail.json').read_text())
        reference=json.loads((folder/'reference_upper.json').read_text())
        return data,tail,F(reference['upper'])-F(14,10000)
    raise ValueError('Unknown bounded molecular case')


def run(names,replay=False):
    rows=[];start=time.monotonic()
    for name in names:
        data,tail,b=load_case(name);path=OUT/f'{name}_global.json';stats=None
        if replay:cert=json.loads(path.read_text())
        else:
            cert,stats=propose(data,tail,b);path.write_text(json.dumps(cert,indent=2)+'\n')
        receipt=check(data,tail,cert);rows.append({'case':name,'discovery':stats,'receipt':receipt,'certificate_bytes':path.stat().st_size})
        print(json.dumps({'case':name,'partition':receipt['partition'],'gap_Ha':receipt['gap_lower_float_Ha'],
            'order':receipt['order'],'actions':receipt['H_actions_per_K_vector']}),flush=True)
    forbidden=[x for x in ('numpy','scipy','cvxpy','pyscf') if x in sys.modules]
    if replay and forbidden:raise AssertionError('Numerical import during exact global replay')
    result={'cases':rows,'wall_seconds':time.monotonic()-start,'numerical_packages_loaded':forbidden}
    (OUT/('global_replay.json' if replay else 'global_discovery.json')).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');parser.add_argument('--cases',nargs='+',default=['h6','fresh_h6_1p6','h8']);a=parser.parse_args();run(a.cases,a.replay)
