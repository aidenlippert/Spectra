"""Combine singlet and all-nonsinglet lower proofs without a spin guess."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time
from experiments.marginal_symbolic import mono,add,scale,product,adj,canonical,decode,hermitian,verified_residual
from research.collective_completion_20260914.spin_replay import setup,alpha_shift
from research.certificate_scaling.spin_twirl import twirl

def spin_squared(m):
    z=add(*(mono(((1,i),(0,i)),F(1 if i%2==0 else -1,2)) for i in range(m)))
    plus=add(*(mono(((1,i),(0,i+1))) for i in range(0,m,2)))
    return add(z,product(z,z),product(canonical(adj(plus)),plus))

def ladder_ideal(m,Y):
    plus=add(*(mono(((1,i),(0,i+1))) for i in range(0,m,2)))
    p=product(plus,Y)
    return add(p,canonical(adj(p)))

def check_sector(data,cert):
    t=time.monotonic();h,hs,delta=setup(data);m,n=data['modes'],data['particles']
    if cert.get('kind')!='spin_sector_sos_v1' or (cert.get('modes'),cert.get('particles'))!=(m,n):raise ValueError('Spin-sector schema')
    mag=cert.get('magnetization');singlet=cert.get('singlet')
    if type(mag) is not int or mag not in (0,1) or type(singlet) is not bool or (mag==1 and singlet):raise ValueError('Spin-sector labels')
    if mag==1 and min(n,m-n)<2:raise ValueError('Empty M_S=1 sector')
    if decode(cert['hamiltonian'],m,4)!=h:raise ValueError('Input Hamiltonian mismatch')
    Y=decode(cert['alpha_multiplier'],m,2)
    if not hermitian(Y) or any(sum(2*c-1 for c,i in w) or sum(2*c-1 for c,i in w if i%2==0) for w in Y):raise ValueError('Spin multiplier conservation')
    if type(cert.get('casimir_multiplier')) is not str:raise ValueError('Exact Casimir coefficient required')
    a=F(cert['casimir_multiplier'])
    if a and not singlet:raise ValueError('Casimir equality requires the singlet sector')
    W=decode(cert.get('spin_ladder_multiplier',[]),m,2)
    if W and not singlet:raise ValueError('Spin ladder equality requires the singlet sector')
    if any(sum(2*c-1 for c,i in w) or sum(2*c-1 for c,i in w if i%2==0)!=-1 for w in W):raise ValueError('Spin ladder multiplier must lower M_S by one')
    Z=add(alpha_shift(m,n),mono((),-mag));expected=add(hs,scale(product(Z,Y),-1),scale(spin_squared(m),-a),scale(ladder_ideal(m,W),-1))
    core=cert['core']
    if core.get('operator_degree')!=3 or (core['modes'],core['particles'])!=(m,n):raise ValueError('Core sector mismatch')
    if decode(core['hamiltonian'],m,4)!=expected:raise ValueError('Auxiliary Hamiltonian identity mismatch')
    residual,rec=verified_residual(core)
    if cert.get('spin_twirl',False):
        if cert.get('spin_twirl') is not True or not singlet:raise ValueError('This spin-average proof requires the singlet sector')
        averaged=twirl(residual);eta=sum(abs(c) for c in averaged.values());bound=F(core['b'])-eta
        rec.update(unaveraged_residual_l1=rec['residual_l1'],residual_l1=str(eta),residual_l1_float=float(eta),
                   lower=str(bound),lower_float=float(bound),residual_terms=len(averaged),
                   spin_twirl=True,positivity='Exact SU(2) average of positive squares; ideals have zero singlet compression')
        residual=averaged
    if 'residual_wedge_witness' in cert:
        from research.collective_completion_20260914.residual_wedge import improve
        rec=improve(core,residual,rec,cert['residual_wedge_witness'])
    rec.update(restricted_to='S=0' if singlet else f'M_S={mag}',original_H_spin_defect_Ha=str(delta),
               lower_for='Exact SU(2) average of input H in the stated sector',replay_seconds=time.monotonic()-t)
    return rec

def check(data,singlet,nonsinglet):
    start=time.monotonic()
    if not singlet.get('singlet') or nonsinglet.get('magnetization')!=1 or nonsinglet.get('singlet'):raise ValueError('Both spin pieces required')
    s=check_sector(data,singlet);t=check_sector(data,nonsinglet);delta=F(s['original_H_spin_defect_Ha'])
    lower=min(F(s['lower']),F(t['lower']))-delta
    return {'kind':'complete_spin_screen_v1','lower':str(lower),'lower_float':float(lower),'singlet':s,'all_nonsinglets':t,
        'valid_on':'Entire fixed-N sector','argument':'Every integer-spin S>=1 multiplet meets M_S=1; S=0 is certified separately; original-H spin defect charged once',
        'many_body_states_enumerated':0,'replay_seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('singlet');p.add_argument('nonsinglet');p.add_argument('output');p.add_argument('--upper');a=p.parse_args()
    read=lambda p:json.loads(Path(p).read_text())
    rec=check(read(a.fixture),read(a.singlet),read(a.nonsinglet))
    if a.upper:
        U=F(read(a.upper)['upper_Ha']);L=F(rec['lower'])
        if U<L:raise ValueError('Inconsistent endpoints')
        rec.update(upper_Ha=str(U),width_Ha=str(U-L),width_mHa=float(1000*(U-L)),target_1p6mHa_met=U-L<=F(1,625))
    if any(k in sys.modules for k in ('numpy','scipy','cvxpy','pyscf','quimb')):raise AssertionError('Numerical acceptance dependency')
    with Path(a.output).open('x') as f:json.dump(rec,f,indent=2)
    print(json.dumps({k:v for k,v in rec.items() if k not in ('singlet','all_nonsinglets','upper_Ha','width_Ha')},indent=2))
