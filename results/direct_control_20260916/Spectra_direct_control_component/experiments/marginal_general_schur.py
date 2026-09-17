"""Exact Schur transfer around a matched reference, with arbitrary real H.

The reference remains the M=10,N=5 matched t=1/5 Hamiltonian. The new H
need only be real, Hermitian, and particle-number conserving (degree <=4).
No pair charges or exchange symmetry of H are assumed. For the signed
symmetric embedding Z, D=Z^T Z and P=Z D^-1 Z^T. A certified reference
complement gap c and eta=sum|coeff(H-H0)| imply QHQ>=c-eta. With
W=(I-P)(H-H0)Z, the exact six-dimensional sufficient test is
    Z^T H Z - bD - W^T W/(c-eta-b) > 0.
Lower construction enumerates reference and coupled configurations; no
large-system scaling claim is made. Floating diagonalization proposes only an upper witness.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import time

from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_schur_transfer import collective_matrix,complement_lower,ldl_pivots
from experiments.marginal_symbolic import add,adj,scale,mono,product,encode,decode,canonical,hermitian,transform
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_symmetry_transfer import hopping_perturbation,integer_upper

ROOT=Path(__file__).resolve().parents[1]


def embedding():
    columns=[{} for _ in range(6)]
    for bits in range(32):
        chosen=[i+5 if bits&(1<<i) else i for i in range(5)]
        inversions=sum(chosen[i]>chosen[j] for i in range(5) for j in range(i+1,5))
        state=sum(1<<i for i in chosen)
        columns[bits.bit_count()][state]=F((-1)**inversions)
    return columns


def apply_columns(h,columns):
    result=[]
    for column in columns:
        output={}
        for word,c in h.items():
            for state,a in column.items():
                target=apply_word(word,state)
                if target:
                    destination,sign=target
                    output[destination]=output.get(destination,F(0))+c*a*sign
        result.append({s:a for s,a in output.items() if a})
    return result


def gram(left,right):
    return [[sum(a*column.get(s,F(0)) for s,a in row.items()) for column in right] for row in left]


def annihilating_moments(h0,columns,max_degree=12):
    """Find and check a common exact recurrence for all coupling columns."""
    return moment_recurrence(columns, lambda current: apply_columns(h0, current), max_degree)


def moment_recurrence(columns, apply, max_degree=12):
    """Exact sparse-vector recurrence for a supplied linear operator action."""
    if type(max_degree) is not int or max_degree < 0:
        raise ValueError('Nonnegative integer recurrence degree budget required')
    echelon=[];powers=[];moments=[];current=columns
    for degree in range(max_degree+1):
        vector={(i,s):a for i,column in enumerate(current) for s,a in column.items() if a}
        relation=[F(0)]*degree+[F(1)]
        for pivot,basis,coefficients in echelon:
            factor=vector.get(pivot,F(0))
            if not factor:continue
            for key,a in basis.items():
                value=vector.get(key,F(0))-factor*a
                if value:vector[key]=value
                else:vector.pop(key,None)
            for i,a in enumerate(coefficients):relation[i]-=factor*a
        if not vector:
            # Independently substitute the relation into the original powers.
            for i in range(len(columns)):
                residual={}
                for coefficient,power in zip(relation,powers+[current]):
                    for s,a in power[i].items():residual[s]=residual.get(s,F(0))+coefficient*a
                if any(residual.values()):raise ValueError('Annihilator substitution failed')
            return {'annihilator':relation,'moments':moments}
        pivot=min(vector);factor=vector[pivot]
        echelon.append((pivot,{key:a/factor for key,a in vector.items()},[a/factor for a in relation]))
        powers.append(current);moments.append(gram(columns,current))
        current=apply(current)
    raise ValueError('No exact annihilator within the configured degree budget')


def resolvent_self_energy(coefficients,moments,z,size=6):
    z=F(z);degree=len(coefficients)-1
    if degree==0:return [[F(0) for _ in range(size)] for _ in range(size)]
    quotient=[F(0)]*degree;quotient[-1]=coefficients[-1]
    for i in range(degree-2,-1,-1):quotient[i]=coefficients[i+1]+z*quotient[i+1]
    value=coefficients[0]+z*quotient[0]
    if not value:raise ValueError('Resolvent parameter is an annihilator root')
    return [[-sum(a*matrix[i][j] for a,matrix in zip(quotient,moments))/value
             for j in range(size)] for i in range(size)]


def norm_bound(delta,method='coefficient_l1'):
    if method=='coefficient_l1':return sum(abs(c) for c in delta.values())
    if method!='hermitian_pairs':raise ValueError('Unknown perturbation norm method')
    # A nondiagonal number-conserving CAR monomial maps between disjoint
    # occupation subspaces. M^2=0, so ||c M + c* M^dagger||<=|c|.
    # A diagonal monomial is a signed occupation projector, norm<=1.
    seen=set();bound=F(0)
    for word,c in delta.items():
        if word in seen:continue
        partner=canonical(adj(mono(word)))
        if len(partner)!=1:raise ValueError('Expected canonical adjoint monomial')
        target,sign=next(iter(partner.items()))
        if delta.get(target)!=c*sign:raise ValueError('Unpaired Hermitian coefficient')
        seen.update((word,target));bound+=abs(c)
    return bound


def prepare(h,norm_method='coefficient_l1',resolvent=False):
    if type(resolvent) is not bool:raise ValueError('Boolean resolvent setting required')
    h=canonical(h)
    if not hermitian(h):raise ValueError('Real Hermitian Hamiltonian required')
    if any(len(word)>4 or sum(1 if creation else -1 for creation,mode in word) for word in h):
        raise ValueError('Particle-number-conserving degree <=4 required')
    if any(not 0<=mode<10 for word in h for creation,mode in word):raise ValueError('Invalid mode')
    h0=hopping_polynomial(10,F(1,5));delta=add(h,scale(h0,-1));z=embedding()
    metric=[F(comb(5,k)) for k in range(6)]
    if gram(z,z)!=[[metric[i] if i==j else F(0) for j in range(6)] for i in range(6)]:
        raise ValueError('Embedding metric mismatch')
    base=apply_columns(h0,z);collective=collective_matrix()
    for j,column in enumerate(base):
        expected={s:a*collective[i][j] for i in range(6) for s,a in z[i].items() if collective[i][j]}
        if column!=expected:raise ValueError('Reference does not preserve signed embedding')
    dz=apply_columns(delta,z);projected=gram(z,dz)
    w=[]
    for j,column in enumerate(dz):
        output=dict(column)
        for i in range(6):
            factor=projected[i][j]/metric[i]
            for s,a in z[i].items():output[s]=output.get(s,F(0))-a*factor
        w.append({s:a for s,a in output.items() if a})
    if any(value for row in gram(z,w) for value in row):raise ValueError('Complement projection failed')
    data={'metric':metric,'projected_h':gram(z,apply_columns(h,z)),
            'embedding':z,'coupling':w,
            'leakage':gram(w,w),'perturbation_norm_bound':norm_bound(delta,norm_method),'norm_method':norm_method,
            'complement_lower':complement_lower(),'embedding_support':sum(len(c) for c in z),
            'leakage_support':len(set().union(*(set(c) for c in w)))}
    if resolvent:data.update(annihilating_moments(h0,w))
    return data


def pivots(data,b):
    b=F(b);q=data['complement_lower']-data['perturbation_norm_bound']-b
    if q<=0:return None
    correction=(resolvent_self_energy(data['annihilator'],data['moments'],b+data['perturbation_norm_bound'])
                if 'annihilator' in data else [[value/q for value in row] for row in data['leakage']])
    matrix=[[data['projected_h'][i][j]-(b*data['metric'][i] if i==j else 0)-correction[i][j]
             for j in range(6)] for i in range(6)]
    return ldl_pivots(matrix)


def lower(data,tolerance=F(1,10**9)):
    tolerance=F(tolerance)
    if not 0<tolerance<1:raise ValueError('Invalid search tolerance')
    hi=min(data['projected_h'][i][i]/data['metric'][i] for i in range(6))
    lo=min(F(3)-data['perturbation_norm_bound'],hi-1);step=F(1)
    for _ in range(64):
        if pivots(data,lo) is not None:break
        lo-=step;step*=2
    else:raise ValueError('Failed to bracket a sufficient lower bound')
    if pivots(data,hi) is not None:raise ValueError('Invalid upper search endpoint')
    while hi-lo>tolerance:
        mid=(hi+lo)/2
        if pivots(data,mid) is not None:lo=mid
        else:hi=mid
    return {'lower':str(lo),'failed_search_endpoint':str(hi),'tolerance':str(tolerance)}


def upper_rayleigh(h,witness):
    states=[s for s in range(1024) if s.bit_count()==5];index={s:i for i,s in enumerate(states)}
    if type(witness) is not dict:raise ValueError('Expected upper witness')
    amps=witness.get('amplitudes')
    if type(amps) is not list or len(amps)!=252 or any(type(a) is not int for a in amps) or not any(amps):
        raise ValueError('Expected integer upper witness in ascending fixed-N basis')
    norm=sum(a*a for a in amps);energy=F(0)
    for word,c in h.items():
        for state,a in zip(states,amps):
            if not a:continue
            target=apply_word(word,state)
            if target:
                destination,sign=target;energy+=c*a*amps[index[destination]]*sign
    return energy/norm


def replay(certificate):
    if certificate.get('kind')!='general_perturbation_schur_v1' or certificate.get('modes')!=10 or certificate.get('particles')!=5:
        raise ValueError('Unsupported certificate')
    h=decode(certificate['hamiltonian'],10,4);data=prepare(h,certificate.get('norm_method','coefficient_l1'),certificate.get('resolvent',False));b=F(certificate['lower'])
    diagonal=pivots(data,b)
    if diagonal is None:raise ValueError('Exact Schur positivity failed')
    upper=upper_rayleigh(h,certificate.get('independent_upper'))
    if upper<b:raise ValueError('Inconsistent interval')
    return {'lower':str(b),'lower_float':float(b),'upper':str(upper),'upper_float':float(upper),
            'width':str(upper-b),'width_float':float(upper-b),'perturbation_norm_bound':str(data['perturbation_norm_bound']),
            'norm_method':data['norm_method'],
            'resolvent_degree':len(data['annihilator'])-1 if 'annihilator' in data else None,
            'reference_complement_lower':str(data['complement_lower']),'schur_dimension':6,
            'embedding_support':data['embedding_support'],'leakage_support':data['leakage_support'],
            'schur_pivots':[str(x) for x in diagonal],'hamiltonian_bound':True,
            'broken_pair_charges':[i for i in range(5) if any(sum((1 if c else -1) for c,m in w if m in (i,i+5)) for w in h)],
            'global_exchange_invariant':transform(h,[(i+5)%10 for i in range(10)])==h,
            'scope':'Real number-conserving degree<=4 H; accuracy is perturbative around a fixed matched reference; no scaling claim'}


def fixture(strength=F(1,1000),interaction=False,modes=10):
    strength=F(strength)
    if strength<=0:raise ValueError('Positive perturbation strength required')
    if type(modes) is not int or modes<4 or modes%2:raise ValueError('Even M>=4 required')
    if interaction and modes<6:raise ValueError('Mixed fixture requires at least three pairs')
    pairs=modes//2
    h=add(hopping_polynomial(modes,F(1,5)),hopping_perturbation(modes,strength)[0])
    # Matched edges plus these cross-pair edges form a connected cycle.
    # Unequal coefficients also remove the fixture's flavor permutations.
    for i in range(pairs):
        p=mono(((1,i),(0,(i+1)%pairs+pairs)))
        h=add(h,scale(add(p,adj(p)),-strength*F(pairs+i,pairs)))
    if interaction:
        transfer=mono(((1,0),(1,1),(0,pairs+2),(0,pairs+1)))
        n0=mono(((1,0),(0,0)));n6=mono(((1,pairs+1),(0,pairs+1)))
        h=add(h,scale(add(transfer,adj(transfer)),strength/2),
              scale(product(n0,n6),strength/7),scale(mono(((1,2),(0,2))),strength/11))
    return canonical(h)


def run(strength,interaction=False,norm_method='coefficient_l1',resolvent=False):
    strength=F(strength);name=('mixed_' if interaction else 'cycle_')+str(strength).replace('/','_')
    out=ROOT/'results/marginal_general_schur'/(name+('_paired' if norm_method=='hermitian_pairs' else '')+('_resolvent' if resolvent else ''))
    if out.exists():raise ValueError('Preserve previous run')
    started=time.monotonic();h=fixture(strength,interaction);data=prepare(h,norm_method,resolvent);search=lower(data)
    certificate={'kind':'general_perturbation_schur_v1','modes':10,'particles':5,
                 'hamiltonian':encode(h),'lower':search['lower'],'norm_method':norm_method,'resolvent':resolvent,'independent_upper':integer_upper(h,10,5)}
    receipt=replay(certificate);receipt['total_seconds']=time.monotonic()-started
    recipe={'strength':str(strength),'interaction':interaction,'norm_method':norm_method,'resolvent':resolvent}
    if resolvent:recipe['annihilator']=[str(x) for x in data['annihilator']]
    out.mkdir(parents=True)
    for key,obj in [('certificate',certificate),('receipt',receipt),('search',search),('fixture',recipe)]:
        (out/(key+'.json')).write_text(json.dumps(obj,indent=2)+'\n')
    print(json.dumps({key:receipt[key] for key in ('lower_float','upper_float','width_float','perturbation_norm_bound','leakage_support','total_seconds')}),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--strength',default='1/1000');parser.add_argument('--interaction',action='store_true');parser.add_argument('--verify',type=Path)
    parser.add_argument('--norm-pairs',action='store_true')
    parser.add_argument('--resolvent',action='store_true')
    args=parser.parse_args()
    if args.verify:print(json.dumps(replay(json.loads(args.verify.read_text())),indent=2))
    else:run(F(args.strength),args.interaction,'hermitian_pairs' if args.norm_pairs else 'coefficient_l1',args.resolvent)
