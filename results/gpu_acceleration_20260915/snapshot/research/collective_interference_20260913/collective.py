"""Implicit bath envelopes with exact retained-space energy certificates.

The retained basis uses unnormalized bright creation operators sum_j b_j†.
Their CAR metric is the bin multiplicity. No bath orbital list or large Fock
sector is made by this module. NumPy is imported only by the proposer.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import hashlib
import json
import time

from experiments.marginal_symbolic import decode,encode,mono,add,scale,product,hermitian
from experiments.marginal_transfer_verify import apply_word
from research.certificate_scaling.commutator_dual_witness import psd


def digest(model):
    return hashlib.sha256(json.dumps(model,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def model(length=2,gap='4',spread='1/2',hybridization='3/5'):
    ns=[mono(((1,i),(0,i))) for i in range(4)]
    h=add(*(scale(ns[i],F(i-2,10)) for i in range(4)))
    for i,j,t in ((0,1,F(-1,2)),(1,2,F(-2,5)),(0,2,F(1,3)),(2,3,F(-3,7)),(0,3,F(-1,5))):
        h=add(h,mono(((1,i),(0,j)),t),mono(((1,j),(0,i)),t))
    for i,j,u in ((0,1,F(1)),(1,2,F(4,5)),(2,3,F(3,5))):h=add(h,scale(product(ns[i],ns[j]),u))
    return {'kind':'two_band_collective_charge_v1','active_modes':4,'active_number_offset':2,
            'bath_length':length,'gap':str(F(gap)),'spread':str(F(spread)),
            'hybridization':str(F(hybridization)),
            'active_hamiltonian':encode(h),'positive_channel':['1','0','1/2','0'],
            'negative_channel':['0','1','0','-1'],
            'charge_u':'1/5','charge_v':'2/7','charge_w':'1/6'}


def parameters(data,bins):
    if data.get('kind')!='two_band_collective_charge_v1' or data.get('active_modes')!=4:
        raise ValueError('Unsupported collective model')
    length=data['bath_length'];offset=data['active_number_offset']
    if type(length) is not int or length<1 or type(offset) is not int or not 0<=offset<=4:
        raise ValueError('Invalid bath multiplicity or fixed-number offset')
    if type(bins) is not int or not 1<=bins<=2:raise ValueError('One or two bins required by retained-mode budget')
    keys=('gap','spread','hybridization','charge_u','charge_v','charge_w')
    if any(not isinstance(data[k],str) for k in keys):raise ValueError('Rational parameter strings required')
    gap,spread,hyb,u,v,w=(F(data[k]) for k in keys)
    if gap<=0 or spread<0:raise ValueError('Positive band gap and nonnegative spread required')
    channels={}
    for sign,key in ((1,'positive_channel'),(-1,'negative_channel')):
        if len(data[key])!=4 or any(not isinstance(x,str) for x in data[key]):raise ValueError('Four rational channel coefficients required')
        channels[sign]=list(map(F,data[key]))
    h=decode(data['active_hamiltonian'],4,4)
    if not hermitian(h) or any(sum(2*c-1 for c,_ in word) for word in h):raise ValueError('Hermitian number-conserving active H required')
    count=length*length;bins=min(bins,count);groups=[]
    for sign in (-1,1):
        for j in range(bins):
            lo=j*count//bins;hi=(j+1)*count//bins-1
            def excitation(index):return gap+(spread*index/(count-1) if count>1 else F(0))
            groups.append({'sign':sign,'count':hi-lo+1,'lo':excitation(lo),'hi':excitation(hi),
                           'index_first':lo,'index_last':hi,'channel':channels[sign]})
    if 4+len(groups)>8:raise ValueError('Retained-mode budget exceeded')
    # On N=count+offset, N_bath-count = -(N_active-offset).
    charge=add(*(mono(((1,i),(0,i))) for i in range(4)),mono((),-offset))
    h=add(h,scale(product(charge,charge),u+v-w))
    filled=-count*(gap+(spread/2 if count>1 else 0))
    return {'h':h,'groups':groups,'count':count,'particles':count+offset,
            'physical_modes':2*count+4,'retained_modes':4+len(groups),
            'filled_energy':filled,'hyb_per_orbital':hyb/length,
            'folded_charge_coefficient':u+v-w}


def states(m,n):
    return [sum(1<<i for i in inds) for inds in combinations(range(m),n)]


def spectator_energy(p,core_particles,endpoint):
    """Minimize all dark occupations by multiplicity, with exact total N."""
    needed=p['particles']-core_particles
    available=sum(g['count']-1 for g in p['groups'])
    if not 0<=needed<=available:return None
    constant=sum(g[endpoint]*(g['count']-1) for g in p['groups'] if g['sign']<0)
    costs=sorted((g['sign']*g[endpoint],g['count']-1) for g in p['groups'])
    value=constant
    for cost,count in costs:
        take=min(needed,count);value+=take*cost;needed-=take
    if needed:raise AssertionError('Spectator count mismatch')
    return value


def core_forms(p,endpoint):
    if endpoint not in ('lo','hi'):raise ValueError('Unknown energy envelope')
    m=p['retained_modes'];gammas=[1]*4+[g['count'] for g in p['groups']]
    h=dict(p['h'])
    for j,g in enumerate(p['groups'],4):
        for i,c in enumerate(g['channel']):
            if c:
                t=c*p['hyb_per_orbital']
                h=add(h,mono(((1,i),(0,j)),t),mono(((1,j),(0,i)),t))
    forms=[];actions=0
    for n in range(m+1):
        dark=spectator_energy(p,n,endpoint)
        if dark is None:continue
        basis=states(m,n);lookup={s:i for i,s in enumerate(basis)}
        metric=[]
        for s in basis:
            value=1
            for i,gamma in enumerate(gammas):
                if s>>i&1:value*=gamma
            metric.append(value)
        matrix=[[F(0) for _ in basis] for _ in basis]
        for col,s in enumerate(basis):
            bath=sum(g[endpoint]*((s>>j&1) if g['sign']>0 else (1-(s>>j&1))) for j,g in enumerate(p['groups'],4))
            matrix[col][col]+=(bath+dark)*metric[col]
            for word,c in h.items():
                target=apply_word(word,s);actions+=1
                if target is None:continue
                t,sign=target;row=lookup[t];factor=1
                for flag,i in word:
                    if not flag:factor*=gammas[i]
                matrix[row][col]+=c*sign*factor*metric[row]
        if any(matrix[i][j]!=matrix[j][i] for i in range(len(basis)) for j in range(len(basis))):
            raise AssertionError('Bright-mode metric Hermiticity failed')
        forms.append({'particles':n,'metric':metric,'matrix':matrix,'spectator_energy':dark})
    return forms,{'retained_states_constructed':sum(len(f['metric']) for f in forms),
                  'largest_matrix':max(len(f['metric']) for f in forms),
                  'retained_word_state_actions':actions,
                  'spectator_bins_processed':len(p['groups'])*len(forms),
                  'physical_fock_space_enumerated':False}


def check_lower(forms,b):
    return [psd([[value-(b*f['metric'][i] if i==j else 0) for j,value in enumerate(row)]
                 for i,row in enumerate(f['matrix'])]) for f in forms]


def rayleigh(forms,witness):
    if type(witness.get('particles')) is not int:raise ValueError('Invalid upper particle sector')
    f=next((f for f in forms if f['particles']==witness['particles']),None)
    if f is None:raise ValueError('Upper sector incompatible with total N')
    v=witness['amplitudes']
    if not isinstance(v,list) or len(v)!=len(f['metric']) or any(type(x) is not int for x in v) or not any(v):
        raise ValueError('Nonzero integer upper coordinates required')
    norm=sum(g*x*x for g,x in zip(f['metric'],v))
    value=sum(v[i]*f['matrix'][i][j]*v[j] for i in range(len(v)) for j in range(len(v)))
    return value/norm


def replay(data,certificate):
    start=time.monotonic()
    if certificate.get('kind')!='collective_bath_envelope_certificate_v1' or certificate.get('model_sha256')!=digest(data):
        raise ValueError('Certificate/model binding failed')
    if not isinstance(certificate['lower_shifted'],str):raise ValueError('Rational lower endpoint required')
    p=parameters(data,certificate['bins']);lower_forms,lower_cost=core_forms(p,'lo');upper_forms,upper_cost=core_forms(p,'hi')
    lower=F(certificate['lower_shifted']);checks=check_lower(lower_forms,lower)
    upper=rayleigh(upper_forms,certificate['upper'])
    low_trial=rayleigh(lower_forms,certificate['lower_trial'])
    active_floor=-sum(abs(c) for c in p['h'].values())
    coupling_norm=2*abs(F(data['hybridization']))*sum(abs(F(x)) for name in ('positive_channel','negative_channel') for x in data[name])
    delta=max(g['hi']-g['lo'] for g in p['groups']);gap=min(g['lo'] for g in p['groups'])
    uniform_upper=low_trial+delta/gap*(low_trial-active_floor+coupling_norm)
    upper=min(upper,uniform_upper)
    if upper<lower:raise ValueError('Inconsistent envelope interval')
    width=upper-lower
    return {'lower':str(lower+p['filled_energy']),'upper':str(upper+p['filled_energy']),
            'lower_shifted':str(lower),'upper_shifted':str(upper),'width':str(width),'width_float':float(width),
            'filled_bath_energy':str(p['filled_energy']),'physical_modes':p['physical_modes'],'particles':p['particles'],
            'retained_modes':p['retained_modes'],'coupling_patterns':len(p['groups']),
            'dark_modes_eliminated':p['physical_modes']-p['retained_modes'],
            'lower_checks':checks,'lower_cost':lower_cost,'upper_cost':upper_cost,
            'maximum_metric_bits':max(g.bit_length() for f in lower_forms for g in f['metric']),
            'uniform_width_ceiling':str(uniform_upper-lower),'lower_envelope_trial_gap':str(low_trial-lower),
            'active_floor':str(active_floor),'hybridization_norm_bound':str(coupling_norm),
            'maximum_bin_width':str(delta),'minimum_excitation_gap':str(gap),
            'replay_seconds':time.monotonic()-start,
            'scope':'Exact operator-envelope interval for declared collective charging/bath family; energy units are abstract, not Hartree.'}


def propose(data,bins=1):
    import numpy as np
    import math
    start=time.monotonic();p=parameters(data,bins)
    low,cost=core_forms(p,'lo');high,_=core_forms(p,'hi')
    def normalized(f):
        g=f['metric']
        return np.array([[float(v/g[i])*math.sqrt(float(F(g[i],g[j]))) for j,v in enumerate(row)] for i,row in enumerate(f['matrix'])])
    minima=[float(np.linalg.eigvalsh(normalized(f))[0]) for f in low]
    lower=F(math.floor(min(minima)*10**10),10**10);attempts=0
    for attempts in range(8):
        try:check_lower(low,lower);break
        except ValueError:lower-=F(10**attempts,10**10)
    else:raise ValueError('Numerical lower proposal failed exact PSD')
    def trial(f):
        _,vectors=np.linalg.eigh(normalized(f))
        v=[float(x)/math.sqrt(g) for x,g in zip(vectors[:,0],f['metric'])]
        rounding=10**10*(math.isqrt(max(f['metric']))+1)
        largest=max(map(abs,v));v=[round(x/largest*rounding) for x in v]
        return {'particles':f['particles'],'amplitudes':v}
    low_trials=[trial(f) for f in low]
    low_value,low_witness=min(((rayleigh(low,w),w) for w in low_trials),key=lambda x:x[0])
    uppers=[]
    for witness in [trial(f) for f in high]+low_trials:
        uppers.append((rayleigh(high,witness),witness))
    value,witness=min(uppers,key=lambda x:x[0])
    certificate={'kind':'collective_bath_envelope_certificate_v1','model_sha256':digest(data),
                 'bins':bins,'lower_shifted':str(lower),'upper':witness,'lower_trial':low_witness}
    return certificate,{'proposal_seconds':time.monotonic()-start,'exact_lower_corrections':attempts,
                        'lower_cost':cost,'numerical_eigensolver_max_dimension':cost['largest_matrix'],
                        'large_sector_enumerated':False}
