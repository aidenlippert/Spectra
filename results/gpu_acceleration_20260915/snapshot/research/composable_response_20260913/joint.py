"""A joint occupation-sector certificate carries the induced interaction intact."""
from fractions import Fraction as F
import argparse
import json
import sys
import time

from experiments.marginal_symbolic import mono,add,product,scale
from research.compact_response_20260913 import program,trial
from research.molecular_collective_20260913.core import extract,tail_replay,digest
from research.positive_response_20260913.coercivity import spatial_inputs,spatial_norm
from research.positive_response_20260913.molecular_diagnostic import spatial_one_body
from research.certificate_scaling.commutator_dual_witness import psd
from research.composable_response_20260913 import spatial_gap,coupling,recursive

OUT=recursive.OUT


def inputs(data,tail,high,x,mu):
    s=data['modes']//2
    if data['modes'] not in (12,16) or data['particles']!=s or type(high) is not int or high not in (2,3):
        raise ValueError('Bounded joint high-orbital sector required')
    if mu<0:raise ValueError('High-occupation multiplier must be nonnegative')
    p=extract(data,tail['center_number']);c,t=spatial_one_body(p);matrices,weights,_=spatial_inputs(data,tail)
    if len(x)!=len(weights):raise ValueError('One tangent per retained density required')
    S=[[t[i][j]+sum(w*v*L[i][j] for w,v,L in zip(weights,x,matrices))-(mu if i==j and i>=s-high else 0)
        for j in range(s)] for i in range(s)]
    base=c+F(tail_replay(data,tail)['lower_operator_shift_Ha'])-sum(w*v*v/2 for w,v in zip(weights,x))+2*mu
    return base,S,p,matrices,weights,t


def check_gap(data,tail,cert):
    start=time.monotonic()
    if cert.get('kind')!='joint_double_occupancy_sector_v1' or cert['fixture_sha256']!=digest(data) or cert['tail_sha256']!=digest(tail):
        raise ValueError('Joint sector input binding failed')
    if any(type(v) is not str for v in cert['tangents']) or any(type(cert[k]) is not str for k in ('occupation_multiplier','chemical_potential','lower_H_joint_Ha')):
        raise ValueError('Exact joint multipliers required')
    base,S,_,_,_,_=inputs(data,tail,cert['last_spatial_orbitals'],list(map(F,cert['tangents'])),F(cert['occupation_multiplier']))
    d=len(S);raw=cert['negative_part_majorant']
    if len(raw)!=d or any(len(row)!=d or any(type(v) is not str for v in row) for row in raw):
        raise ValueError('Wrong joint orbital majorant shape')
    Y=[list(map(F,row)) for row in raw];tau=F(cert['chemical_potential'])
    psd(Y);psd([[S[i][j]+Y[i][j]-(tau if i==j else 0) for j in range(d)] for i in range(d)])
    gamma=base+data['particles']*tau-2*sum(Y[i][i] for i in range(d))
    if F(cert['lower_H_joint_Ha'])>gamma:raise ValueError('Joint endpoint overstated')
    return {'lower_H_joint_Ha':cert['lower_H_joint_Ha'],'proved_lower_H_joint_Ha':str(gamma),
        'largest_matrix_dimension':d*(d+1)//2,'new_orbital_matrix_dimension':d,
        'tail_coefficient_matrix_dimension':d*(d+1)//2,'new_orbital_PSD_checks':2,'many_body_states_enumerated':0,
        'many_body_matrix_entries':0,'replay_seconds':time.monotonic()-start}


def propose_gap(data,tail,high=2):
    import numpy as np
    from scipy.optimize import minimize
    start=time.monotonic();matrices,weights,_=spatial_inputs(data,tail);s=data['modes']//2;n=data['particles']
    base,S,_,_,_,_=inputs(data,tail,high,[F(0)]*len(weights),F(0))
    w=np.array(list(map(float,weights)));Ls=np.array(matrices,dtype=float);T=np.array(S,dtype=float);P=np.diag([0.]*(s-high)+[1.]*high)
    def objective(y):
        x=y[:-1];mu=y[-1];e,U=np.linalg.eigh(T+np.einsum('k,kij->ij',w*x,Ls)-mu*P);occ=U[:,:n//2]
        value=float(base)-np.sum(w*x*x)/2+2*mu+2*np.sum(e[:n//2])
        gradient=np.r_[w*(2*np.einsum('ia,kij,ja->k',occ,Ls,occ)-x),2-2*np.sum(occ[-high:,:]**2)]
        return -value,-gradient
    result=minimize(objective,np.zeros(len(w)+1),jac=True,method='L-BFGS-B',bounds=[(None,None)]*len(w)+[(0,None)],
        options={'maxiter':500,'ftol':1e-13,'gtol':1e-9})
    x=[F(round(float(v)*10**8),10**8) for v in result.x[:-1]];mu=F(round(float(result.x[-1])*10**8),10**8)
    base,S,_,_,_,_=inputs(data,tail,high,x,mu);e,U=np.linalg.eigh(np.array(S,dtype=float));k=n//2
    tau=F(round(float((e[k-1]+e[k])/2)*10**8),10**8)
    Yf=(U*np.maximum(float(tau)-e,0.))@U.T+1e-8*np.eye(s)
    Y=[[F(round(float((Yf[i,j]+Yf[j,i])/2)*10**10),10**10) for j in range(s)] for i in range(s)]
    gamma=base+n*tau-2*sum(Y[i][i] for i in range(s))
    cert={'kind':'joint_double_occupancy_sector_v1','fixture_sha256':digest(data),'tail_sha256':digest(tail),
        'last_spatial_orbitals':high,'tangents':list(map(str,x)),'occupation_multiplier':str(mu),
        'chemical_potential':str(tau),'negative_part_majorant':[[str(v) for v in row] for row in Y],
        'lower_H_joint_Ha':str(program.floor_grid(gamma,10**10))}
    receipt=check_gap(data,tail,cert)
    return cert,{'construction_seconds':time.monotonic()-start,'optimizer_iterations':int(result.nit),
        'optimizer_success':bool(result.success),'optimizer_message':str(result.message),'acceptance':receipt}


def second_endpoints(data,tail,first,first_gap_cert,joint_cert):
    first_receipt=program.check(data,tail,first);b=F(first['target_Ha']);eta=F(first_receipt['exact_program_residual_penalty_Ha'])
    gamma=F(check_gap(data,tail,joint_cert)['lower_H_joint_Ha']);joint_delta=gamma-b-eta
    delta1=F(spatial_gap.check(data,tail,first_gap_cert)['lower_H_QQ_Ha'])-b
    if delta1<=0:raise ValueError('Positive first denominator bound required')
    s=data['modes']//2;j=s-2;n=data['particles']-2;rem=[i for i in range(s) if i!=j]
    p=extract(data,tail['center_number']);constant,t=spatial_one_body(p);matrices,weights,_=spatial_inputs(data,tail)
    upper=constant+2*t[j][j]+program.one_body_endpoint([[t[i][k] for k in rem] for i in rem],n,True)
    for w,L in zip(weights,matrices):
        A=[[L[i][k] for k in rem] for i in rem];vv=sum(L[i][j]**2 for i in rem)
        upper+=w*((n*spatial_norm(A)+2*abs(L[j][j]))**2/2+vv)
    M2=upper+F(tail_replay(data,tail)['upper_operator_shift_Ha'])-b-eta
    q=product(mono(((1,2*j),(0,2*j))),mono(((1,2*j+1),(0,2*j+1))))
    bare_coupling=product(product(q,p['h']),add(mono(()),scale(q,-1)))
    bare_g=sum(abs(v) for v in bare_coupling.values());cross=coupling.prove(data)
    g2=bare_g+coupling.sqrt_up(F(cross['coupling_norm_squared_Ha2']))*F(first['coupling_norm_Ha'])/delta1
    return {'delta':joint_delta,'M':M2,'g':g2,'gamma_joint':gamma,'delta1':delta1,'first_eta':eta,
        'bare_coupling_terms':len(bare_coupling),'local_coupling':cross}


def outer_program(bounds):
    delta=program.floor_grid(bounds['delta']);M=program.ceil_grid(bounds['M']);g=program.ceil_grid(bounds['g']);budget=F(1,10**6)
    if not 0<delta<M:return {'status':'no_positive_second_gap','delta_Ha':str(delta)}
    z=(M+delta)/(M-delta)
    for k in range(1,257):
        eta=g*g/(delta*program.chebyshev(k,z)**2)
        if eta<=budget:
            return {'status':'certified_response','delta_Ha':str(delta),'M_Ha':str(M),'coupling_norm_Ha':str(g),
                'order':k,'residual_budget_Ha':str(budget),'residual_penalty_Ha':str(eta)}
    return {'status':'order_cap_exhausted','delta_Ha':str(delta),'M_Ha':str(M),'coupling_norm_Ha':str(g),'order_cap':256,
        'penalty_at_cap_Ha':str(eta)}


def load_case(name):
    if name=='fresh_h6_1p6':
        source=program.ROOT/'results/response_consistency_20260913/fresh_h6_1p6'
        data=json.loads((source/'fixture.json').read_text());tail=json.loads((source/'tail.json').read_text())
        first=json.loads((OUT/'transfer_program.json').read_text())
    else:
        if name not in ('h6','h8'):raise ValueError('Only the three declared molecular fixtures are allowed')
        source=program.ROOT/'results/molecular_collective_20260913/campaign'/name
        rank={'h6':10,'h8':14}[name]
        data=json.loads((source/'fixture.json').read_text());tail=json.loads((source/f'rank_{rank}/tail.json').read_text())
        first=json.loads((program.OUT/f'{name}_program.json').read_text())
    # The lower-bound component needs b from the response program, not the
    # inherited enumerated upper amplitudes. Their discovery is charged separately.
    return data,tail,None,first


def check(data,tail,first,cert):
    start=time.monotonic()
    if cert.get('kind')!='composed_joint_sector_response_v1' or cert['fixture_sha256']!=digest(data) or cert['first_response_sha256']!=digest(first):
        raise ValueError('Composed response binding failed')
    if cert['joint_sector']['last_spatial_orbitals']!=2 or cert['first_sector']['orbital']!=data['modes']//2-1:
        raise ValueError('The composed response requires the declared two-level partition')
    bounds=second_endpoints(data,tail,first,cert['first_sector'],cert['joint_sector']);expected=outer_program(bounds)
    if cert['second_response']!=expected:raise ValueError('Composed response scalar program did not reproduce')
    status=expected['status'];k1=first['order'];k2=expected.get('order')
    return {'status':status,'target_Ha':first['target_Ha'],'joint_gap_before_first_penalty_Ha':str(bounds['gamma_joint']-F(first['target_Ha'])),
        'second_sector_gap_Ha':str(bounds['delta']),'second_sector_gap_float_Ha':float(bounds['delta']),
        'second_M_float_Ha':float(bounds['M']),'second_coupling_norm_float_Ha':float(bounds['g']),
        'first_order':k1,'second_order':k2,'base_H_actions_per_second_response_right_hand_side':(k2-1)*(2*k1+1) if k2 else None,
        'base_H_actions_per_twice_retained_operator_vector':(2*k2+1)*(2*k1+1) if k2 else None,
        'second_residual_penalty_Ha':expected.get('residual_penalty_Ha'),'many_body_states_enumerated':0,'many_body_matrix_entries':0,
        'largest_orbital_matrix':data['modes']//2,
        'tail_coefficient_matrix_dimension':(data['modes']//2)*(data['modes']//2+1)//2,
        'local_transition_labels':7,'coupling_work':bounds['local_coupling'],
        'replay_seconds':time.monotonic()-start,'terminal_retained_positivity_certified':False,
        'scope':'The second response and its gap are certified without determinant enumeration. Final retained positivity and physical costs of nested operator actions remain separate.'}


def run(replay=False):
    start=time.monotonic();rows=[]
    for name in ('h6','fresh_h6_1p6','h8'):
        data,tail,reference,first=load_case(name);path=OUT/f'{name}_joint_response.json';discovery=None
        if replay:cert=json.loads(path.read_text())
        else:
            t=time.monotonic();gc,gs=propose_gap(data,tail);fc,fs=spatial_gap.propose(data,tail,data['modes']//2-1)
            bounds=second_endpoints(data,tail,first,fc,gc)
            cert={'kind':'composed_joint_sector_response_v1','fixture_sha256':digest(data),'first_response_sha256':digest(first),
                'joint_sector':gc,'first_sector':fc,'second_response':outer_program(bounds)}
            path.write_text(json.dumps(cert,indent=2)+'\n');discovery={'joint_gap':gs,'first_gap':fs,'construction_seconds':time.monotonic()-t}
        receipt=check(data,tail,first,cert);rows.append({'case':name,'receipt':receipt,'discovery':discovery,'certificate_bytes':path.stat().st_size})
        print(json.dumps({'case':name,'status':receipt['status'],'second_gap_Ha':receipt['second_sector_gap_float_Ha'],'second_order':receipt['second_order']}),flush=True)
    result={'cases':rows,'wall_seconds':time.monotonic()-start,'numerical_packages_loaded':[x for x in ('numpy','scipy','cvxpy','pyscf') if x in sys.modules]}
    if replay and result['numerical_packages_loaded']:raise AssertionError('Numerical import during composed response replay')
    (OUT/('joint_replay.json' if replay else 'joint_discovery.json')).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');run(parser.parse_args().replay)
