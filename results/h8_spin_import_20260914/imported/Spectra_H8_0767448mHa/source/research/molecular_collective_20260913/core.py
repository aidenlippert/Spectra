"""Hamiltonian-derived density patterns and a collective tail certificate.

All acceptance uses the standard library. Numerical packages propose factors,
PSD shifts, supporting hyperplanes, and an exact-tested dual density only.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import json
import time

from experiments.marginal_symbolic import canonical,decode,encode,mono,add,scale,product,hermitian
from experiments.marginal_general_schur import norm_bound
from experiments.marginal_transfer_verify import apply_word
from research.certificate_scaling.commutator_dual_witness import psd


def digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def matmul(a,b):
    return [[sum((x*b[k][j] for k,x in enumerate(row)),F(0)) for j in range(len(b[0]))] for row in a]


def transpose(a):return [list(row) for row in zip(*a)]


def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]


def coefficient(h,word):
    image=canonical(mono(tuple(word)))
    if len(image)!=1:raise ValueError('Unique canonical monomial required')
    w,sign=next(iter(image.items()))
    return h.get(w,F(0))/sign


def density(p,q):
    return add(*(mono(((1,2*p+s),(0,2*q+s))) for s in range(2)))


def square_form(matrix,basis):
    out={}
    for i,row in enumerate(matrix):
        for j,c in enumerate(row):
            if c:
                for w,a in product(basis[i],basis[j]).items():out[w]=out.get(w,F(0))+c*a/2
    return {w:c for w,c in out.items() if c}


def extract(data,center_number=True):
    start=time.monotonic();m=data['modes'];n=data['particles']
    if type(m) is not int or m%2 or not 4<=m<=20 or type(n) is not int or not 0<=n<=m:
        raise ValueError('Two through ten spatial orbitals and a valid fixed particle number required')
    if type(center_number) is not bool:raise ValueError('Boolean number-centering flag required')
    s=m//2;h=decode(data['hamiltonian'],m,4)
    if not hermitian(h) or any(sum(2*c-1 for c,_ in w) for w in h):raise ValueError('Hermitian number-conserving H required')
    pairs=list(combinations_with_replacement(range(s),2))
    basis=[density(p,q) if p==q else add(density(p,q),density(q,p)) for p,q in pairs]
    metric=[F(1) if p==q else F(1,2) for p,q in pairs]
    def eri(p,q,r,t):
        return coefficient(h,((1,2*p),(1,2*r+1),(0,2*t+1),(0,2*q)))
    matrix=[]
    for p,q in pairs:
        row=[]
        for r,t in pairs:
            row.append(sum((eri(a,b,c,d) for a,b in ((p,q),(q,p)) for c,d in ((r,t),(t,r))),F(0))/8+
                       sum((eri(c,d,a,b) for a,b in ((p,q),(q,p)) for c,d in ((r,t),(t,r))),F(0))/8)
        matrix.append(row)
    whole=square_form(matrix,basis)
    remainder=add(h,scale(whole,-1))
    one={w:c for w,c in remainder.items() if len(w)<=2}
    mismatch=add(remainder,scale(one,-1))
    eta=norm_bound(mismatch,'hermitian_pairs')
    casimir=F(n*(s+2))-F(n*n,2)
    if center_number:
        # Q=P Q + (N/s)u. Fold every term containing the conserved N into
        # the one-body part on the requested sector before factor discovery.
        u=[F(i==j) for i,j in pairs];cu=[sum(c*v for c,v in zip(row,u)) for row in matrix]
        total=sum(x*y for x,y in zip(u,cu))
        matrix=[[c-u[i]*cu[j]/s-cu[i]*u[j]/s+u[i]*u[j]*total/(s*s)
                 for j,c in enumerate(row)] for i,row in enumerate(matrix)]
        one=add(one,*(scale(q,F(n,s)*(c-v*total/s)) for q,c,v in zip(basis,cu,u)),
                mono((),F(n*n,2*s*s)*total))
        casimir-=F(n*n,s)
    return {'h':h,'modes':m,'particles':n,'spatial':s,'pairs':pairs,'basis':basis,'metric':metric,
            'matrix':matrix,'one':one,'mismatch':mismatch,'eta':eta,'casimir':casimir,
            'center_number':center_number,'extract_seconds':time.monotonic()-start}


def parse_factors(items,d):
    if not isinstance(items,list) or len(items)>d:raise ValueError('Factor count exceeds feature budget')
    result=[]
    for item in items:
        if not isinstance(item['weight'],str) or F(item['weight'])<0:raise ValueError('Nonnegative rational factor weight required')
        v=item['vector']
        if not isinstance(v,list) or len(v)!=d or any(not isinstance(c,str) for c in v):raise ValueError('Rational factor coordinates required')
        v=list(map(F,v))
        if not any(v):raise ValueError('Zero factor vector')
        result.append((F(item['weight']),v))
    return result


def tail_matrix(p,factors):
    d=len(p['basis'])
    return [[p['matrix'][i][j]-sum((w*v[i]*v[j] for w,v in factors),F(0)) for j in range(d)] for i in range(d)]


def tail_replay(data,cert):
    start=time.monotonic()
    if cert.get('kind')!='molecular_density_tail_v1' or cert.get('fixture_sha256')!=digest(data):raise ValueError('Fixture binding failed')
    p=extract(data,cert['center_number']);d=len(p['basis']);factors=parse_factors(cert['factors'],d)
    if p['center_number'] and any(sum(c for c,(i,j) in zip(v,p['pairs']) if i==j) for _,v in factors):
        raise ValueError('Centered factors must have exactly zero orbital trace')
    if any(not isinstance(cert[k],str) for k in ('alpha','beta')):raise ValueError('Rational tail endpoints required')
    alpha=F(cert['alpha']);beta=F(cert['beta'])
    if alpha>0 or beta<0:raise ValueError('Tail envelope must straddle zero')
    residual=tail_matrix(p,factors)
    lower=psd([[c-(alpha*p['metric'][i] if i==j else 0) for j,c in enumerate(row)] for i,row in enumerate(residual)])
    upper=psd([[(beta*p['metric'][i] if i==j else 0)-c for j,c in enumerate(row)] for i,row in enumerate(residual)])
    lo=alpha*p['casimir']/2-p['eta'];hi=beta*p['casimir']/2+p['eta']
    return {'lower_operator_shift_Ha':str(lo),'upper_operator_shift_Ha':str(hi),
            'tail_interval_width_Ha':str(hi-lo),'tail_interval_width_mHa':float(1000*(hi-lo)),
            'factors':len(factors),'factor_nonzeros':sum(sum(bool(x) for x in v) for _,v in factors),
            'coefficient_dimension':d,'modes':p['modes'],'particles':p['particles'],
            'number_centered':p['center_number'],
            'mismatch_paired_norm_Ha':str(p['eta']),'mismatch_terms':len(p['mismatch']),
            'casimir_bound':str(p['casimir']),'lower_psd':lower,'upper_psd':upper,
            'replay_seconds':time.monotonic()-start,'many_body_states_enumerated':0,
            'scope':'Operator interval H_retained + lo I <= H <= H_retained + hi I; no retained-system energy solve implied.'}


def propose_tail(data,p,rank):
    import math
    import numpy as np
    start=time.monotonic();d=len(p['basis'])
    if type(rank) is not int or not 0<=rank<=d:raise ValueError('Invalid retained rank')
    g=np.sqrt(np.array(p['metric'],dtype=float));matrix=np.array(p['matrix'],dtype=float)/g[:,None]/g[None,:]
    values,vectors=np.linalg.eigh(matrix);factors=[]
    for index in range(d-1,d-rank-1,-1):
        weight=F(max(0,round(float(values[index])*10**12)),10**12)
        vector=[F(round(float(x)*10**10),10**10) for x in g*vectors[:,index]]
        if p['center_number']:
            mean=sum(c for c,(i,j) in zip(vector,p['pairs']) if i==j)/p['spatial']
            vector=[c-mean if i==j else c for c,(i,j) in zip(vector,p['pairs'])]
        if weight:factors.append({'weight':str(weight),'vector':[str(x) for x in vector]})
    residual=tail_matrix(p,parse_factors(factors,d))
    vals=np.linalg.eigvalsh(np.array(residual,dtype=float)/g[:,None]/g[None,:])
    alpha=F(min(0,math.floor(float(vals[0])*10**10))-1,10**10)
    beta=F(max(0,math.ceil(float(vals[-1])*10**10))+1,10**10)
    cert={'kind':'molecular_density_tail_v1','fixture_sha256':digest(data),'center_number':p['center_number'],
          'factors':factors,'alpha':str(alpha),'beta':str(beta)}
    for retry in range(7):
        try:
            psd([[c-(alpha*p['metric'][i] if i==j else 0) for j,c in enumerate(row)] for i,row in enumerate(residual)])
            psd([[(beta*p['metric'][i] if i==j else 0)-c for j,c in enumerate(row)] for i,row in enumerate(residual)])
            break
        except ValueError:
            alpha-=F(10**retry,10**10);beta+=F(10**retry,10**10)
    else:raise ValueError('Tail spectral proposal failed exact verification')
    cert.update(alpha=str(alpha),beta=str(beta))
    return cert,{'proposal_seconds':time.monotonic()-start,'eigenvalues_descending':values[::-1].tolist(),
                 'exact_shift_retries':retry,'largest_numeric_matrix':d,'physical_sector_enumerated':False}


def retained_polynomial(p,cert):
    factors=parse_factors(cert['factors'],len(p['basis']));h=dict(p['one'])
    for weight,vector in factors:
        q=add(*(scale(op,c) for op,c in zip(p['basis'],vector) if c))
        h=add(h,scale(product(q,q),weight/2))
    return h


def one_matrix(poly,m):
    out=[[F(0) for _ in range(m)] for _ in range(m)]
    for w,c in poly.items():
        if not w:continue
        if len(w)!=2 or (w[0][0],w[1][0])!=(1,0):raise ValueError('One-body polynomial required')
        out[w[0][1]][w[1][1]]+=c
    if out!=transpose(out):raise ValueError('Symmetric one-body matrix required')
    return out


def factor_operators(p,cert):
    return [(w,add(*(scale(q,c) for q,c in zip(p['basis'],v) if c)))
            for w,v in parse_factors(cert['factors'],len(p['basis']))]


def rational_matrix(raw,m):
    if not isinstance(raw,list) or len(raw)!=m or any(not isinstance(row,list) or len(row)!=m for row in raw):
        raise ValueError('Invalid exact matrix shape')
    if any(not isinstance(c,str) for row in raw for c in row):raise ValueError('Rational matrix strings required')
    return [list(map(F,row)) for row in raw]


def pairing(a,b):return sum((a[i][j]*b[j][i] for i in range(len(a)) for j in range(len(a))),F(0))


def support_replay(data,tail,cert):
    start=time.monotonic();receipt=tail_replay(data,tail);p=extract(data,tail['center_number']);m=p['modes'];n=p['particles']
    if cert.get('kind')!='density_square_support_v1' or cert.get('tail_sha256')!=digest(tail):raise ValueError('Supporting-bound binding failed')
    factors=factor_operators(p,tail)
    if len(cert['centers'])!=len(factors) or any(not isinstance(x,str) for x in cert['centers']):raise ValueError('Rational centers required')
    centers=list(map(F,cert['centers']));linear=dict(p['one']);constant=F(0)
    for (weight,q),center in zip(factors,centers):
        linear=add(linear,scale(q,weight*center));constant-=weight*center*center/2
    k=one_matrix(linear,m);mu=F(cert['chemical_potential']);shift=F(cert['psd_shift'])
    if shift<0:raise ValueError('Nonnegative PSD correction required')
    v=rational_matrix(cert['negative_factor'],m);x=matmul(v,transpose(v))
    for i in range(m):x[i][i]+=shift
    psd([[k[i][j]+x[i][j]-(mu if i==j else 0) for j in range(m)] for i in range(m)])
    retained_lower=linear.get((),F(0))+constant+n*mu-sum(x[i][i] for i in range(m))
    # Feasible one-body density gives a ceiling for every choice of centers.
    gamma=rational_matrix(cert['dual_density'],m)
    if sum(gamma[i][i] for i in range(m))!=n:raise ValueError('Dual trace differs from particle number')
    psd(gamma);psd([[F(i==j)-gamma[i][j] for j in range(m)] for i in range(m)])
    ceiling=p['one'].get((),F(0))+pairing(one_matrix(p['one'],m),gamma)
    for weight,q in factors:ceiling+=weight*pairing(one_matrix(q,m),gamma)**2/2
    if ceiling<retained_lower:raise ValueError('Supporting primal/dual ordering failed')
    tail_lower=F(receipt['lower_operator_shift_Ha']);lower=retained_lower+tail_lower
    hf=(1<<n)-1;upper=F(0)
    for w,c in p['h'].items():
        target=apply_word(w,hf)
        if target is not None and target[0]==hf:upper+=c*target[1]
    return {'certified_lower_Ha':str(lower),'hf_upper_Ha':str(upper),'interval_width_Ha':str(upper-lower),
            'retained_lower_Ha':str(retained_lower),'support_family_ceiling_Ha':str(ceiling+tail_lower),
            'support_primal_dual_gap_Ha':str(ceiling-retained_lower),'factors':len(factors),
            'many_body_states_enumerated':0,'upper_determinants_evaluated':1,'largest_support_psd':m,
            'replay_seconds':time.monotonic()-start,
            'scope':'Ground-energy lower from linear supports of retained squares plus tail lower; exact dual ceiling only for this supporting-bound grammar.'}


def propose_support(data,p,tail):
    import numpy as np
    import cvxpy as cp
    start=time.monotonic();m=p['modes'];n=p['particles'];factors=factor_operators(p,tail)
    h=np.array(one_matrix(p['one'],m),float);operators=[np.array(one_matrix(q,m),float) for _,q in factors]
    gamma=cp.Variable((m,m),symmetric=True)
    objective=cp.trace(h@gamma)
    for (weight,_),q in zip(factors,operators):objective+=float(weight)/2*cp.square(cp.trace(q@gamma))
    problem=cp.Problem(cp.Minimize(objective),[gamma>>0,np.eye(m)-gamma>>0,cp.trace(gamma)==n])
    problem.solve(solver='CLARABEL',max_iter=200,tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
    if gamma.value is None or not np.all(np.isfinite(gamma.value)):raise ValueError('Supporting-bound proposal failed')
    centers=[F(round(float(np.trace(q@gamma.value))*10**8),10**8) for q in operators]
    linear=dict(p['one'])
    for (weight,q),center in zip(factors,centers):linear=add(linear,scale(q,weight*center))
    k=one_matrix(linear,m);ev,u=np.linalg.eigh(np.array(k,float))
    chemical=ev[0]-1 if n==0 else ev[-1]+1 if n==m else (ev[n-1]+ev[n])/2
    mu=F(round(float(chemical)*10**10),10**10)
    vnum=u@np.diag(np.sqrt(np.maximum(float(mu)-ev,0)))
    v=[[F(round(float(c)*10**10),10**10) for c in row] for row in vnum];x=matmul(v,transpose(v))
    shift=F(1,10**9)
    for _ in range(8):
        try:psd([[k[i][j]+x[i][j]+((shift-mu) if i==j else 0) for j in range(m)] for i in range(m)]);break
        except ValueError:shift*=10
    else:raise ValueError('One-body PSD rounding failed')
    raw=[[F(round(float(c)*10**8),10**8) for c in row] for row in (gamma.value+gamma.value.T)/2]
    trace=sum(raw[i][i] for i in range(m));raw[-1][-1]+=n-trace
    mix=F(1,10**7)
    for _ in range(7):
        exact=[[(1-mix)*raw[i][j]+(mix*F(n,m) if i==j else 0) for j in range(m)] for i in range(m)]
        try:psd(exact);psd([[F(i==j)-exact[i][j] for j in range(m)] for i in range(m)]);break
        except ValueError:mix*=10
    else:raise ValueError('Dual density rounding failed')
    cert={'kind':'density_square_support_v1','tail_sha256':digest(tail),'centers':[str(x) for x in centers],
          'chemical_potential':str(mu),'psd_shift':str(shift),
          'negative_factor':[[str(c) for c in row] for row in v],
          'dual_density':[[str(c) for c in row] for row in exact]}
    return cert,{'proposal_seconds':time.monotonic()-start,'solver_status':problem.status,
                 'solver_solve_seconds':problem.solver_stats.solve_time,'dual_interior_mix':str(mix),
                 'physical_sector_enumerated':False,'largest_numeric_one_body_matrix':m}


def rank_replay(data,cert):
    start=time.monotonic()
    if cert.get('kind')!='density_tail_rank_obstruction_v1' or cert.get('fixture_sha256')!=digest(data):raise ValueError('Rank witness binding failed')
    p=extract(data,cert['center_number']);d=len(p['basis']);budget=F(cert['budget_Ha'])
    if budget<=2*p['eta'] or p['casimir']<=0:raise ValueError('Positive nontrivial tail budget required')
    z=cert['subspace'];r=cert['minimum_factors']
    if type(r) is not int or not 1<=r<=d or len(z)!=d or any(len(row)!=r for row in z):raise ValueError('Invalid rank subspace')
    if any(not isinstance(c,str) for row in z for c in row):raise ValueError('Exact rank subspace required')
    z=[list(map(F,row)) for row in z];threshold=2*(budget-2*p['eta'])/p['casimir']
    shifted=[[c-(threshold*p['metric'][i] if i==j else 0) for j,c in enumerate(row)] for i,row in enumerate(p['matrix'])]
    check=psd(matmul(transpose(z),matmul(shifted,z)))
    if check['rank']!=r:raise ValueError('Strict positive rank witness required')
    return {'minimum_factors':r,'budget_Ha':str(budget),'threshold':str(threshold),
            'positive_minor':check,'replay_seconds':time.monotonic()-start,'many_body_states_enumerated':0,
            'scope':'Minimum rank for the coefficient Loewner envelope and stated Casimir width budget; not a no-go for other operator bounds or nonlinear pattern descriptions.'}


def propose_rank(data,p,minimum,budget=F(16,10000)):
    import numpy as np
    g=np.sqrt(np.array(p['metric'],float));values,vectors=np.linalg.eigh(np.array(p['matrix'],float)/g[:,None]/g[None,:])
    z=vectors[:,-minimum:]/g[:,None]
    return {'kind':'density_tail_rank_obstruction_v1','fixture_sha256':digest(data),
            'center_number':p['center_number'],'budget_Ha':str(budget),'minimum_factors':minimum,
            'subspace':[[str(F(round(float(c)*10**8),10**8)) for c in row] for row in z]}
