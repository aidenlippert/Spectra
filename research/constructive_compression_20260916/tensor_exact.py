"""Integer MPS contractions and exact recurrence residuals, no configurations.

SVD, QR and floating point are absent. A tensor network is accepted only as
explicit rational data; its cost follows the network bonds and sparsity.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import prod,isqrt,lcm


def state(cert):
    den=cert['denominator'];q=cert['bond_charges'];layers=cert['tensors'];m=cert['modes']
    if type(den)!=int or den<1 or len(q)!=m+1 or len(layers)!=m:raise ValueError('MPS data')
    if q[0]!=[[0,0]] or q[-1]!=[cert['spin_counts']]:raise ValueError('Sector boundary')
    for i,layer in enumerate(layers):
        seen=set()
        for a,s,b,v in layer:
            if any(type(x)!=int for x in (a,s,b,v)) or s not in (0,1) or not v:raise ValueError('Integer entry')
            if not (0<=a<len(q[i]) and 0<=b<len(q[i+1])) or (a,s,b) in seen:raise ValueError('Index')
            seen.add((a,s,b));qq=list(q[i][a]);qq[i%2]+=s
            if qq!=q[i+1][b]:raise ValueError('Charge flow')
    return {'widths':list(map(len,q)),'denominators':[den]*m,'layers':layers}


def apply_mpo(a,op):
    if len(a['layers'])!=len(op['layers']):raise ValueError('MPO length')
    out=[];dens=[]
    for i,(edges,ops) in enumerate(zip(a['layers'],op['layers'])):
        dd=lcm(*(F(e[3]).denominator for e in ops));new=defaultdict(int)
        for ol,orr,mat,c in ops:
            co=int(F(c)*dd)
            for l,s,r,v in edges:
                for p in (0,1):
                    if mat[2*p+s]:new[ol*a['widths'][i]+l,p,orr*a['widths'][i+1]+r]+=co*mat[2*p+s]*v
        out.append([[l,s,r,v] for (l,s,r),v in sorted(new.items()) if v]);dens.append(dd*a['denominators'][i])
    return {'widths':[x*y for x,y in zip(a['widths'],op['widths'])],'denominators':dens,'layers':out}


def inner(a,b,stats=None):
    if len(a['layers'])!=len(b['layers']):raise ValueError('Inner length')
    e={(0,0):1};work=0;peak=1
    for ea,eb in zip(a['layers'],b['layers']):
        ra=[defaultdict(list),defaultdict(list)];rb=[defaultdict(list),defaultdict(list)]
        for l,s,r,v in ea:ra[s][l].append((r,v))
        for l,s,r,v in eb:rb[s][l].append((r,v))
        out=defaultdict(int)
        for s in (0,1):
            mid=defaultdict(int)
            for (l,k),ev in e.items():
                for r,v in ra[s].get(l,()):mid[r,k]+=v*ev;work+=1
            for (r,k),ev in mid.items():
                if ev:
                    for t,v in rb[s].get(k,()):out[r,t]+=ev*v;work+=1
        e={k:v for k,v in out.items() if v};peak=max(peak,len(e))
    if stats is not None:
        stats['integer_multiplications']=stats.get('integer_multiplications',0)+work
        stats['maximum_transfer_entries']=max(stats.get('maximum_transfer_entries',0),peak)
        stats['contractions']=stats.get('contractions',0)+1
    return F(e.get((0,0),0),prod(a['denominators'])*prod(b['denominators']))


def squared_norm_combination(terms,stats=None):
    total=F(0)
    for i,(c,a) in enumerate(terms):
        for j in range(i+1):
            d,b=terms[j];total+=c*d*inner(a,b,stats)*(1 if i==j else 2)
    if total<0:raise ArithmeticError('Negative exact squared norm')
    return total


def sqrt_upper(x,bits=40):
    x=F(x)
    if x<0 or type(bits)!=int or bits<0:raise ValueError('Square root premise')
    numerator=x.numerator<<(2*bits);q=isqrt(numerator//x.denominator)
    if q*q*x.denominator<numerator:q+=1
    return F(q,1<<bits)


def recurrence_residual(v,previous,older,op,b,ell,stats=None):
    b,ell=F(b),F(ell)
    if b<=ell:raise ValueError('Polynomial scaling')
    hp=apply_mpo(previous,op);alpha=(b+ell)/(b-ell);beta=-2/(b-ell)
    terms=[(F(1),v),(-alpha if older is None else -2*alpha,previous),
           (-beta if older is None else -2*beta,hp)]
    if older is not None:terms.append((F(1),older))
    n2=squared_norm_combination(terms,stats)
    return n2,sqrt_upper(n2)
