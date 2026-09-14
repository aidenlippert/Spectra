"""Exact local-mixture upper limits for a fixed joint-projector relaxation.

An accepted result caps attainable lower certificates, not physical ground
energy. Pure integer vectors with positive rational mixture weights make local
PSD automatic; exact expectations enforce the dual constraints.
"""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _actions, _exact, _sector
from experiments.marginal_projector_extendibility import _projector_vector
from experiments.marginal_charged_projectors import charged_vectors


def _profiles(x):
    a,b,p,q,d,e=x
    return [a,b,10-a-b,10-a-b,b,a],[p,q,5-2*p-2*q,q,p],[d,e,F(5,2)-2*d-2*e,e,d]


def _physical_vector(source):
    if type(source) is not dict or not 1<=len(source)<=400:
        raise ValueError('Bounded local mixture vector required')
    result={}
    for key,a in source.items():
        if type(key) not in (str,int):raise ValueError('Integer determinant required')
        s=int(key)
        if s in result or not 0<=s<4096 or type(a) is not int or abs(a)>10**9:
            raise ValueError('Bounded integer local amplitudes required')
        if a:result[s]=a
    norm=sum(a*a for a in result.values())
    if not norm:raise ValueError('Nonzero mixture source required')
    return result,norm


def _weight(value):
    if type(value) not in (str,int,F) or len(str(value))>4096:
        raise ValueError('Bounded rational mixture weight required')
    w=F(value)
    if not 0<=w<=1:raise ValueError('Nonnegative normalized mixture weights required')
    return w


def replay(c):
    if type(c) is not dict or c.get('kind')!='joint_projector_family_limit_v1':
        raise ValueError('Unsupported joint family limit')
    half,hn=_projector_vector(c.get('half_vector'),6)
    if {_sector(s,6) for s,a in half.items() if a}!={(3,3)}:
        raise ValueError('Half projector must occupy spin sector (3,3)')
    charged,cn=charged_vectors(c.get('charged_vector'))
    ratio=_exact(c.get('ratio'),10**9)
    th=_exact(c.get('theta_half'),10**12);tj=_exact(c.get('theta_joint'),10**12)
    if ratio<=0 or not 0<=th<=1 or not 0<=tj<=max(1,ratio):
        raise ValueError('Invalid fixed-projector ceiling or ratio')
    mixture=c.get('mixture')
    if type(mixture) is not list or not 1<=len(mixture)<=16:
        raise ValueError('One through16 physical mixture sources required')
    states=[]
    for item in mixture:
        if type(item) is not dict:raise ValueError('Explicit physical mixture required')
        v,n=_physical_vector(item.get('vector'));states.append((_weight(item.get('weight')),v,n))
    if sum(w for w,_,_ in states)!=1:raise ValueError('Mixture trace must equal one exactly')
    def fidelity(source,norm):
        return sum((w*F(sum(a*source.get(s,0) for s,a in v.items())**2,n*norm)
                    for w,v,n in states),F(0))
    ph=fidelity(half,hn);q=ph+ratio*sum((fidelity(v,cn) for v in charged),F(0))
    if ph>th or q>tj:raise ValueError('Mixture violates a fixed-projector fidelity ceiling')
    energies=[]
    for index in range(7):
        x=[F(int(index==i+1)) for i in range(6)];u,t,d=_profiles(x)
        actions=_actions(6,F(10,3),1,u,t,F(1,2),d)
        energies.append(sum((w*sum((F(a)*b*v.get(target,0) for s,a in v.items()
                                      for target,b in actions[s].items()),F(0))/n
                             for w,v,n in states),F(0)))
    gradients=[e-energies[0] for e in energies[1:]]
    if any(gradients):raise ValueError('Profile expectations must cancel exactly')
    upper=energies[0]/5
    if 'proposed_periodic_family_upper' in c and F(c['proposed_periodic_family_upper'])!=upper:
        raise ValueError('Proposed family limit disagrees with physical expectation')
    return {'accepted':True,'periodic_family_upper':str(upper),'half_expectation':str(ph),
            'joint_expectation':str(q),'theta_half':str(th),'theta_joint':str(tj),
            'profile_gradients':list(map(str,gradients)),'mixture_trace':'1','mixture_sources':len(states),
            'scope':'Upper limit on attainable periodic lower density in the fixed-source six-site U4,t1,V1/2 joint-projector family with nonnegative alpha,beta. All reflected mean-correct profiles, even signed coefficients, are covered because six exact profile derivatives vanish. This local PSD mixture need not extend globally and is not a physical ground-energy upper bound. Supplied projector ceilings define the relaxation; their all-state validity is certified separately.'}
