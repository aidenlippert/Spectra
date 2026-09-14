"""Exact fixed-N positive-amplitude bounds for open spinless fermion chains.

All accepting arithmetic is rational and standard-library only. The method
requires nearest-neighbor nonpositive hopping in the occupation basis. It
certifies the supplied finite model, not generic molecular accuracy.
"""
from fractions import Fraction as F
from itertools import product
import argparse
import json
import time
from pathlib import Path


def rational(x):
    if type(x) not in (int,str,F):
        raise ValueError('Exact rational input required')
    return F(x)


def validate(model, amplitude):
    m,n=model['modes'],model['particles']
    if type(m) is not int or m<2 or type(n) is not int or not 0<=n<=m:
        raise ValueError('Invalid fermionic sector')
    if m*(n+1)*8>5000000:
        raise ValueError('Declared DP work exceeds this runner budget')
    arrays=[]
    for obj,key,length,positive in ((model,'hopping',m-1,False),(model,'interaction',m-1,None),
                                   (model,'fields',m,None),(amplitude,'sites',m,True),
                                   (amplitude,'bonds',m-1,True)):
        raw=obj[key]
        if not isinstance(raw,list) or len(raw)!=length:
            raise ValueError('Wrong array length: '+key)
        a=[rational(x) for x in raw]
        if positive is True and any(x<=0 for x in a):raise ValueError('Amplitude weights must be strictly positive')
        if positive is False and any(x<0 for x in a):raise ValueError('Nonstoquastic hopping refused')
        arrays.append(a)
    return m,n,*arrays,rational(model.get('offset',0))


def factor_tables(model, amplitude):
    m,n,t,v,h,y,x,offset=validate(model,amplitude)
    ending=[[] for _ in range(m)]
    for i in range(m):ending[i].append((1,[F(0),h[i]]))
    entries=2*m
    for i in range(m-1):
        start=max(0,i-1);end=min(m-1,i+2);width=end-start+1;values=[]
        for code in range(1<<width):
            bit=lambda j:(code>>(end-j))&1
            a,b=bit(i),bit(i+1);e=v[i]*a*b
            if a!=b:
                ratio=(y[i+1]/y[i])**(a-b)
                if i>0:ratio*=x[i-1]**(bit(i-1)*(b-a))
                if i+2<m:ratio*=x[i+1]**(bit(i+2)*(a-b))
                e-=t[i]*ratio
            values.append(e)
        entries+=len(values);ending[end].append((width,values))
    return ending,entries


def bounds(model, amplitude):
    start=time.monotonic();m,n,t,v,h,y,x,offset=validate(model,amplitude)
    ending,entries=factor_tables(model,amplitude)
    # key=(number so far,last 3 bits); values=(minimum,minimum witness,Z,sum(w*E)).
    states={(0,0):(offset,0,F(1),offset)}
    peak_states=1;transitions=0;factor_lookups=0;peak_bits=1
    for i in range(m):
        nxt={};remaining=m-i-1
        for (count,tail),(minimum,witness,Z,S) in states.items():
            for bit in (0,1):
                number=count+bit
                if number>n or number+remaining<n:continue
                transitions+=1;full=(tail<<1)|bit
                local=sum((table[full&((1<<width)-1)] for width,table in ending[i]),F(0))
                factor_lookups+=len(ending[i])
                weight=y[i]**(2*bit)
                if i:weight*=x[i-1]**(2*(tail&1)*bit)
                z=Z*weight;s=(S+Z*local)*weight
                lo=minimum+local;mask=(witness<<1)|bit;key=(number,full&7)
                if key in nxt:
                    oldlo,oldmask,oldz,olds=nxt[key]
                    if (oldlo,oldmask)<(lo,mask):lo,mask=oldlo,oldmask
                    z+=oldz;s+=olds
                nxt[key]=(lo,mask,z,s)
                peak_bits=max(peak_bits,*(max(a.numerator.bit_length(),a.denominator.bit_length()) for a in (lo,z,s)))
        states=nxt;peak_states=max(peak_states,len(states))
    lower,mask=min((a[0],a[1]) for a in states.values())
    Z=sum((a[2] for a in states.values()),F(0));S=sum((a[3] for a in states.values()),F(0))
    if Z<=0:raise AssertionError('Nonpositive normalization')
    upper=S/Z
    if upper<lower:raise AssertionError('Inconsistent interval')
    return {'method':'fixed_N_positive_Jastrow_local_energy_DP_v1','lower':str(lower),'upper':str(upper),
            'width':str(upper-lower),'lower_float':float(lower),'upper_float':float(upper),
            'width_float':float(upper-lower),'partition':str(Z),
            'minimum_configuration':format(mask,f'0{m}b'),'modes':m,'particles':n,
            'dp_transitions':transitions,'peak_dp_states':peak_states,'factor_table_entries':entries,
            'factor_lookups':factor_lookups,'peak_rational_bits':peak_bits,'wall_seconds':time.monotonic()-start,
            'many_body_states_enumerated':0,'scope':'Exact rational bounds for specified open nearest-neighbor chain; no generic accuracy guarantee.'}


def local_energy(model, amplitude, configuration):
    """Evaluate the declared local factor sum, useful for replay witnesses."""
    m,n,*_=validate(model,amplitude)
    if len(configuration)!=m or any(type(b) is not int or b not in (0,1) for b in configuration) or sum(configuration)!=n:
        raise ValueError('Invalid fixed-N configuration')
    ending,_=factor_tables(model,amplitude);tail=0;e=rational(model.get('offset',0))
    for i,bit in enumerate(configuration):
        tail=(tail<<1)|bit
        e+=sum((table[tail&((1<<width)-1)] for width,table in ending[i]),F(0))
    return e


def replay(payload):
    actual=bounds(payload['model'],payload['amplitude'])
    if 'claim' in payload:
        for key in ('lower','upper','width','partition','minimum_configuration'):
            if payload['claim'][key]!=actual[key]:raise ValueError('False stored claim: '+key)
    return actual


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--witness',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args();r=replay(json.loads(a.witness.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
