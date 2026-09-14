"""Exact Schur certificate with automatically retained coupling states.

Expand the signed symmetric reference space by the exact H0-Krylov closure
of Q delta Z. Treat the actual Hamiltonian on that expanded space exactly;
bound only its remaining leakage. A fixed reference complement gap still
supplies the final bound. No orthonormal or floating basis is trusted.
"""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,lcm
import json
import time

from experiments.marginal_general_schur import prepare as base_prepare,apply_columns,gram,fixture,upper_rayleigh,annihilating_moments,resolvent_self_energy
from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_symbolic import encode,decode
from experiments.marginal_symmetry_transfer import integer_upper

ROOT=Path(__file__).resolve().parents[1]


def reduced(column,echelon):
    result=dict(column)
    for pivot,basis in echelon:
        factor=result.get(pivot,F(0))
        if not factor:continue
        for s,a in basis.items():
            value=result.get(s,F(0))-factor*a
            if value:result[s]=value
            else:result.pop(s,None)
    return result


def closure(h0,columns,max_dimension=64,orthogonal=False):
    if orthogonal:
        basis=[];norms=[];queue=list(columns);index=0
        def remove(column):
            vector=dict(column)
            for b,norm in zip(basis,norms):
                factor=sum(a*vector.get(s,F(0)) for s,a in b.items())/norm
                if not factor:continue
                for s,a in b.items():
                    value=vector.get(s,F(0))-factor*a
                    if value:vector[s]=value
                    else:vector.pop(s,None)
            return vector
        while index<len(queue):
            vector=remove(queue[index]);index+=1
            if not vector:continue
            if len(basis)>=max_dimension:raise ValueError('Krylov dimension budget exhausted')
            factor=max(abs(a) for a in vector.values());vector={s:a/factor for s,a in vector.items()}
            basis.append(vector);norms.append(sum(a*a for a in vector.values()))
            queue.extend(apply_columns(h0,[vector]))
        if any(remove(column) for column in apply_columns(h0,basis)):raise ValueError('Orthogonal closure check failed')
        return basis
    echelon=[];queue=list(columns);index=0
    while index<len(queue):
        vector=reduced(queue[index],echelon);index+=1
        if not vector:continue
        if len(echelon)>=max_dimension:raise ValueError('Krylov dimension budget exhausted')
        pivot=min(vector);factor=vector[pivot];vector={s:a/factor for s,a in vector.items()}
        echelon.append((pivot,vector));queue.extend(apply_columns(h0,[vector]))
    basis=[v for _,v in echelon]
    if any(reduced(column,echelon) for column in apply_columns(h0,basis)):
        raise ValueError('Reference closure verification failed')
    return basis


def solve_positive(matrix,rhs):
    n=len(matrix)
    if not n or len(rhs)!=n or any(len(row)!=n for row in matrix):raise ValueError('Invalid linear system')
    a=[[F(x) for x in row] for row in matrix];b=[[F(x) for x in row] for row in rhs]
    if any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):raise ValueError('Nonsymmetric metric')
    width=len(b[0])
    if any(len(row)!=width for row in b):raise ValueError('Inconsistent right-hand sides')
    if not any(a[i][j] for i in range(n) for j in range(n) if i!=j):
        if any(a[i][i]<=0 for i in range(n)):raise ValueError('Nonpositive metric pivot')
        return [[value/a[i][i] for value in b[i]] for i in range(n)]
    for i in range(n):
        d=a[i][i]-sum(a[i][k]**2*a[k][k] for k in range(i))
        if d<=0:raise ValueError('Nonpositive metric pivot')
        for j in range(i+1,n):a[j][i]=(a[j][i]-sum(a[j][k]*a[k][k]*a[i][k] for k in range(i)))/d
        a[i][i]=d
    y=[[F(0)]*width for _ in range(n)]
    for i in range(n):
        for j in range(width):y[i][j]=b[i][j]-sum(a[i][k]*y[k][j] for k in range(i))
    x=[[F(0)]*width for _ in range(n)]
    for i in range(n-1,-1,-1):
        for j in range(width):x[i][j]=y[i][j]/a[i][i]-sum(a[k][i]*x[k][j] for k in range(i+1,n))
    if any(sum(matrix[i][k]*x[k][j] for k in range(n))!=rhs[i][j] for i in range(n) for j in range(width)):
        raise ValueError('Exact metric solve residual')
    return x


def projected_action(h,u):
    metric=gram(u,u);hu=apply_columns(h,u);projected=gram(u,hu)
    coefficients=solve_positive(metric,projected);leak=[]
    for j,column in enumerate(hu):
        output=dict(column)
        for i,basis in enumerate(u):
            for s,a in basis.items():output[s]=output.get(s,F(0))-a*coefficients[i][j]
        leak.append({s:a for s,a in output.items() if a})
    if any(a for row in gram(u,leak) for a in row):raise ValueError('Leakage projection failed')
    return metric,projected,leak


def prepare(h,resolvent=False,rounds=1,target_coefficients=None,orthogonal=False):
    if type(resolvent) is not bool:raise ValueError('Boolean resolvent setting required')
    if type(orthogonal) is not bool:raise ValueError('Boolean orthogonal setting required')
    if type(rounds) is not int or not 1<=rounds<=2:raise ValueError('One or two enrichment rounds required')
    if target_coefficients is not None and rounds!=2:raise ValueError('Target requires a second round')
    started=time.monotonic();base=base_prepare(h,'hermitian_pairs');h0=hopping_polynomial(10,F(1,5))
    z=base['embedding'];u=list(z);seeds=base['coupling'];dimensions=[len(u)]
    for iteration in range(rounds):
        if iteration==1 and target_coefficients is not None:
            if type(target_coefficients) is not list or len(target_coefficients)!=len(seeds) or any(type(a) is not int for a in target_coefficients) or not any(target_coefficients):
                raise ValueError('Expected nonzero integer targeting coefficients')
            selected={}
            for a,column in zip(target_coefficients,seeds):
                for s,value in column.items():selected[s]=selected.get(s,F(0))+a*value
            seeds=[{s:a for s,a in selected.items() if a}]
        r=closure(h0,seeds,max_dimension=64-len(u),orthogonal=orthogonal)
        if any(a for row in gram(u,r) for a in row):raise ValueError('New Krylov space is not orthogonal to retained space')
        u.extend(r);dimensions.append(len(u));metric,projected,leak=projected_action(h,u);seeds=leak
    if any(leak[:6]):raise ValueError('Initial coupling was not retained exactly')
    if orthogonal and any(metric[i][j] for i in range(len(u)) for j in range(len(u)) if i!=j):
        raise ValueError('Exact orthogonal metric check failed')
    data={'metric':metric,'projected_h':projected,'leakage':gram(leak,leak),
            'basis_columns':u,'hamiltonian':h,'orthogonal':orthogonal,
            'norm_bound':base['perturbation_norm_bound'],'complement_lower':base['complement_lower'],
            'retained_dimension':len(u),'coupling_dimension':len(u)-6,'stage_dimensions':dimensions,
            'retained_support':len(set().union(*(set(c) for c in u))),
            'leakage_support':len(set().union(*(set(c) for c in leak)))}
    if resolvent:data.update(annihilating_moments(h0,leak))
    data['build_seconds']=time.monotonic()-started
    return data


def pivots(data,b):
    b=F(b);q=data['complement_lower']-data['norm_bound']-b
    if q<=0:return None
    n=data['retained_dimension']
    correction=(resolvent_self_energy(data['annihilator'],data['moments'],b+data['norm_bound'],size=n)
                if 'annihilator' in data else [[value/q for value in row] for row in data['leakage']])
    matrix=[[data['projected_h'][i][j]-b*data['metric'][i][j]-correction[i][j]
             for j in range(n)] for i in range(n)]
    return ldl_pivots(matrix)


def lower(data,tolerance=F(1,10**9)):
    tolerance=F(tolerance)
    if not 0<tolerance<1:raise ValueError('Invalid search tolerance')
    n=data['retained_dimension'];hi=min(data['projected_h'][i][i]/data['metric'][i][i] for i in range(n))
    lo=min(F(3)-data['norm_bound'],hi-1);step=F(1)
    for _ in range(64):
        if pivots(data,lo) is not None:break
        lo-=step;step*=2
    else:raise ValueError('Failed lower bracketing')
    if pivots(data,hi) is not None:raise ValueError('Invalid failed endpoint')
    started=time.monotonic();checks=0
    while hi-lo>tolerance:
        middle=(hi+lo)/2;checks+=1
        if pivots(data,middle) is not None:lo=middle
        else:hi=middle
    return {'lower':str(lo),'failed_search_endpoint':str(hi),'tolerance':str(tolerance),
            'search_seconds':time.monotonic()-started,'positivity_checks':checks}


def physical_coordinates(data):
    import numpy as np
    states=[s for s in range(1024) if s.bit_count()==5]
    columns=np.array([[float(column.get(s,0)) for column in data['basis_columns']] for s in states])
    columns/=np.linalg.norm(columns,axis=0)
    q,_=np.linalg.qr(columns)
    actions=apply_columns(data['hamiltonian'],[{t:F(1)} for t in states])
    h=np.array([[float(column.get(s,0)) for column in actions] for s in states])
    # The full fixed-sector matrix is used only by the numerical proposer.
    return states,q,h


def suggest_target(data):
    import numpy as np
    from scipy.linalg import eigh
    h=np.array(data['projected_h'],dtype=float);g=np.array(data['metric'],dtype=float)
    if data['orthogonal']:
        values,vectors=eigh(h,g,subset_by_index=[0,0]);vector=vectors[:,0];largest=max(abs(vector))
        coefficients=[int(round(float(a/largest)*10000)) for a in vector]
        proposal={'selector':'orthogonal_coordinates','integer_scale':10000}
    else:
        states,q,full=physical_coordinates(data);values,vectors=eigh(q.T@full@q,subset_by_index=[0,0])
        vector=q@vectors[:,0];largest=max(abs(vector))
        physical={s:F(round(float(a/largest)*10**6)) for s,a in zip(states,vector) if round(float(a/largest)*10**6)}
        coordinates=solve_positive(data['metric'],gram(data['basis_columns'],[physical]))
        denominator=lcm(*(row[0].denominator for row in coordinates))
        coefficients=[int(row[0]*denominator) for row in coordinates];divisor=gcd(*coefficients)
        coefficients=[a//divisor for a in coefficients]
        proposal={'selector':'rounded_physical_projection','integer_scale':10**6,
                  'physical_amplitudes':[int(physical.get(s,0)) for s in states],
                  'maximum_coordinate_bits':max(abs(a).bit_length() for a in coefficients)}
    exact_energy=sum(coefficients[i]*data['projected_h'][i][j]*coefficients[j] for i in range(len(g)) for j in range(len(g)))
    exact_norm=sum(coefficients[i]*data['metric'][i][j]*coefficients[j] for i in range(len(g)) for j in range(len(g)))
    return coefficients,dict(proposal,numerical_ritz_energy=float(values[0]),
                         rounded_rayleigh=str(exact_energy/exact_norm),rounded_rayleigh_float=float(exact_energy/exact_norm))


def proposed_lower(data):
    """Use floating search only to propose b; exact positivity accepts it."""
    import numpy as np
    from scipy.linalg import eigvalsh
    import math
    started=time.monotonic();g=np.array(data['metric'],dtype=float)
    scaling=np.sqrt(np.diag(g));outer=np.outer(scaling,scaling);g=g/outer
    h=np.array(data['projected_h'],dtype=float)/outer;k=np.array(data['leakage'],dtype=float)/outer
    c=float(data['complement_lower']);eta=float(data['norm_bound'])
    if 'annihilator' in data:
        f=np.array(data['annihilator'],dtype=float)
        moments=[np.array(matrix,dtype=float)/outer for matrix in data['moments']]
    if 'basis_columns' in data:
        states,orthonormal,full=physical_coordinates(data);action=full@orthonormal
        h=orthonormal.T@action;leak=action-orthonormal@h;k=leak.T@leak;g=np.eye(len(h))
        if 'annihilator' in data:
            actions=apply_columns(hopping_polynomial(10,F(1,5)),[{t:F(1)} for t in states])
            h0=np.array([[float(column.get(s,0)) for column in actions] for s in states])
            moments=[];power=leak
            for _ in range(len(f)-1):moments.append(leak.T@power);power=h0@power
    def positive(b):
        q=c-eta-b
        if q<=0:return False
        correction=k/q
        if 'annihilator' in data:
            degree=len(f)-1
            if degree:
                quotient=np.zeros(degree);quotient[-1]=f[-1]
                for i in range(degree-2,-1,-1):quotient[i]=f[i+1]+(b+eta)*quotient[i+1]
                value=f[0]+(b+eta)*quotient[0]
                correction=-sum(a*m for a,m in zip(quotient,moments))/value
            else:correction=np.zeros_like(h)
        return eigvalsh(h-b*g-correction,g,subset_by_index=[0,0])[0]>0
    hi=float(min(data['projected_h'][i][i]/data['metric'][i][i] for i in range(len(g))))
    lo=min(3-eta,hi-1);step=1.
    for _ in range(64):
        if positive(lo):break
        lo-=step;step*=2
    else:raise ValueError('Numerical proposal failed to bracket')
    for _ in range(64):
        middle=(lo+hi)/2
        if positive(middle):lo=middle
        else:hi=middle
    if not math.isfinite(lo):raise ValueError('Nonfinite proposed lower')
    for check in range(8):
        margin=10.**(check-9);candidate=F(math.floor((lo-margin)*10**10),10**10)
        if pivots(data,candidate) is not None:
            return {'lower':str(candidate),'numerical_search_lower':lo,'numerical_search_upper':hi,
                    'accepted_margin':margin,'exact_positivity_checks':check+1,
                    'search_seconds':time.monotonic()-started,'acceptance':'exact rational positivity'}
    raise ValueError('No exact certificate accepted from numerical search')


def replay(certificate):
    if certificate.get('kind')!='enlarged_coupling_schur_v1' or certificate.get('modes')!=10 or certificate.get('particles')!=5:
        raise ValueError('Unsupported enlarged certificate')
    h=decode(certificate['hamiltonian'],10,4);data=prepare(h,certificate.get('resolvent',False),certificate.get('rounds',1),certificate.get('target_coefficients'),certificate.get('orthogonal',False));b=F(certificate['lower']);diagonal=pivots(data,b)
    if diagonal is None:raise ValueError('Exact enlarged Schur positivity failed')
    upper=upper_rayleigh(h,certificate.get('independent_upper'))
    if upper<b:raise ValueError('Inconsistent interval')
    return {'lower':str(b),'lower_float':float(b),'upper':str(upper),'upper_float':float(upper),
            'width':str(upper-b),'width_float':float(upper-b),'retained_dimension':data['retained_dimension'],
            'coupling_dimension':data['coupling_dimension'],'retained_support':data['retained_support'],
            'stage_dimensions':data['stage_dimensions'],
            'leakage_support':data['leakage_support'],'perturbation_norm_bound':str(data['norm_bound']),
            'resolvent_degree':len(data['annihilator'])-1 if 'annihilator' in data else None,
            'schur_pivots':[str(x) for x in diagonal],'hamiltonian_bound':True,
            'scope':'Fixed matched reference, exact retained coupling space, bounded residual leakage; no general scaling claim'}


def run(strength,interaction=False,resolvent=False,rounds=1,targeted=False,orthogonal=False):
    if targeted:rounds=2
    strength=F(strength);name=('mixed_' if interaction else 'cycle_')+str(strength).replace('/','_')
    out=ROOT/'results/marginal_enlarged_schur'/(name+('_resolvent' if resolvent else '')+('_rounds2' if rounds==2 else '')+('_targeted' if targeted else '')+('_orthogonal' if orthogonal else ('_physical' if targeted else '')))
    if out.exists():raise ValueError('Preserve previous run')
    out.mkdir(parents=True);started=time.monotonic();h=fixture(strength,interaction);coefficients=None
    if targeted:
        coefficients,proposal=suggest_target(prepare(h,orthogonal=orthogonal))
        (out/'target_proposal.json').write_text(json.dumps(dict(proposal,coefficients=coefficients),indent=2)+'\n')
    data=prepare(h,resolvent,rounds,coefficients,orthogonal)
    structure={k:data[k] for k in ('retained_dimension','coupling_dimension','retained_support','leakage_support','build_seconds')}
    (out/'structure.json').write_text(json.dumps(structure,indent=2)+'\n');print(json.dumps(structure),flush=True)
    search=proposed_lower(data) if targeted else lower(data)
    certificate={'kind':'enlarged_coupling_schur_v1','modes':10,'particles':5,'hamiltonian':encode(h),
                 'lower':search['lower'],'resolvent':resolvent,'rounds':rounds,'orthogonal':orthogonal,'independent_upper':integer_upper(h,10,5)}
    if targeted:certificate['target_coefficients']=coefficients
    receipt=replay(certificate);receipt['total_seconds']=time.monotonic()-started
    for key,obj in [('certificate',certificate),('receipt',receipt),('search',search),('fixture',{'strength':str(strength),'interaction':interaction,'resolvent':resolvent,'rounds':rounds,'targeted':targeted,'orthogonal':orthogonal})]:
        (out/(key+'.json')).write_text(json.dumps(obj,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('lower_float','upper_float','width_float','retained_dimension','total_seconds')}),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--strength',default='1/100');parser.add_argument('--interaction',action='store_true');parser.add_argument('--verify',type=Path)
    parser.add_argument('--resolvent',action='store_true')
    parser.add_argument('--rounds',type=int,choices=(1,2),default=1)
    parser.add_argument('--targeted',action='store_true')
    parser.add_argument('--orthogonal',action='store_true')
    args=parser.parse_args()
    if args.verify:print(json.dumps(replay(json.loads(args.verify.read_text())),indent=2))
    else:run(F(args.strength),args.interaction,args.resolvent,args.rounds,args.targeted,args.orthogonal)
