"""Outward binary64 interval contractions of the unchanged rational MPS.

Every rounded operation is enclosed by adjacent floats. Fast math and
parallel reassociation are disabled. Finite IEEE-754 binary64 arithmetic with
gradual underflow is an explicit platform requirement, checked at startup.
Ambiguity never accepts an endpoint. The integer checker is the fallback.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys
import time
import numpy as np
from numba import njit
from experiments.marginal_symbolic import decode,hermitian,add,scale,mono
from research.correlated_pair_20260913.mps_exact import State,I,local_word,check as exact_check
from research.molecular_collective_20260913.core import digest

@njit(cache=True,fastmath=False)
def plus(al,ah,bl,bh):
    if bl==0. and bh==0.:return al,ah
    if al==0. and ah==0.:return bl,bh
    return np.nextafter(al+bl,-np.inf),np.nextafter(ah+bh,np.inf)

@njit(cache=True,fastmath=False)
def times(al,ah,bl,bh):
    if (al==0. and ah==0.) or (bl==0. and bh==0.):return 0.,0.
    a=al*bl;b=al*bh;c=ah*bl;d=ah*bh
    return np.nextafter(min(a,b,c,d),-np.inf),np.nextafter(max(a,b,c,d),np.inf)

@njit(cache=True,fastmath=False)
def transfer(lo,hi,edges,values,O,right):
    left=lo.shape[0];outl=np.zeros((right,right));outh=np.zeros((right,right))
    for s in range(2):
        if O[2*s]==0 and O[2*s+1]==0:continue
        ml=np.zeros((right,left));mh=np.zeros((right,left))
        for k in range(len(edges)):
            a,p,b=edges[k]
            if p!=s:continue
            vl,vh=values[k]
            for j in range(left):
                if lo[a,j]==0. and hi[a,j]==0.:continue
                pl,ph=times(vl,vh,lo[a,j],hi[a,j])
                ml[b,j],mh[b,j]=plus(ml[b,j],mh[b,j],pl,ph)
        for k in range(len(edges)):
            a,p,b=edges[k];ov=O[2*s+p]
            if ov==0:continue
            vl,vh=values[k]
            for j in range(right):
                if ml[j,a]==0. and mh[j,a]==0.:continue
                pl,ph=times(ml[j,a],mh[j,a],vl,vh)
                if ov==-1:pl,ph=-ph,-pl
                outl[j,b],outh[j,b]=plus(outl[j,b],outh[j,b],pl,ph)
    return outl,outh

def rational_interval(q):
    q=F(q);v=float(q)
    if not math.isfinite(v):raise ValueError('Nonfinite coefficient conversion')
    exact=F(v)
    return (math.nextafter(v,-math.inf) if exact>q else v,
            math.nextafter(v,math.inf) if exact<q else v)

def platform_gate():
    if sys.float_info.radix!=2 or sys.float_info.mant_dig!=53:raise ValueError('Binary64 required')
    tiny=math.ldexp(1.,-1022)
    l,h=times(tiny,tiny,0.5,0.5)
    if not F(l)<=F(tiny)/2<=F(h) or l==0:raise ValueError('Gradual underflow required')
    l,h=plus(1.,1.,2.**-54,2.**-54)
    if not F(l)<=1+F(1,2**54)<=F(h):raise ValueError('Outward arithmetic unavailable')

class EnclosedState(State):
    def contract_factors(self,factors):
        # Reuse all structural gates, replacing only the old exact norm call.
        # No caller accepts this sentinel: rigorous norm positivity is checked
        # unconditionally in __init__ below before this object is usable.
        return 1

    def __init__(self,data,cert):
        super().__init__(data,cert)
        self.edges=[np.array([[a,s,b] for a,s,b,v in site],dtype=np.int64).reshape(-1,3) for site in cert['tensors']]
        self.values=[np.array([rational_interval(F(v,self.den)) for a,s,b,v in site],dtype=float).reshape(-1,2) for site in cert['tensors']]
        self.dims=list(map(len,cert['bond_charges']));self.steps=0
        self.norm=self.contract([I]*self.m)
        if not all(math.isfinite(v) for v in self.norm) or self.norm[0]<=0:raise ValueError('Norm positivity not established')

    def step_interval(self,lo,hi,site,op):
        self.steps+=1
        l,h=transfer(lo,hi,self.edges[site],self.values[site],np.array(op,dtype=np.int64),self.dims[site+1])
        if not np.isfinite(l).all() or not np.isfinite(h).all():raise ValueError('Nonfinite enclosure')
        return l,h

    def contract(self,factors):
        lo=hi=np.ones((1,1))
        for site,op in enumerate(factors):lo,hi=self.step_interval(lo,hi,site,op)
        return float(lo[0,0]),float(hi[0,0])

    def polynomial(self,poly):
        trie={}
        for word,c in poly.items():
            node=trie
            for op in local_word(word,self.m):node=node.setdefault(op,{})
            node[None]=node.get(None,F(0))+c
        def walk(node,lo,hi,site):
            if site==self.m:
                return times(float(lo[0,0]),float(hi[0,0]),*rational_interval(node.get(None,F(0))))
            total=(0.,0.)
            for op,child in node.items():
                l,h=self.step_interval(lo,hi,site,op)
                total=plus(*total,*walk(child,l,h,site+1))
            return total
        return walk(trie,np.ones((1,1)),np.ones((1,1)),0) if trie else (0.,0.)

def check(data,cert,endpoint,fallback=False):
    start=time.monotonic();platform_gate();u=F(endpoint)
    state=EnclosedState(data,cert);h=decode(data['hamiltonian'],data['modes'],4)
    if not hermitian(h) or any(sum(2*c-1 for c,i in w) for w in h):raise ValueError('Hermitian number-conserving H required')
    gap=state.polynomial(add(mono((),u),scale(h,F(-1))))
    if not all(math.isfinite(v) for v in gap):raise ValueError('Nonfinite endpoint enclosure')
    rec={'method':'outward_binary64_charge_MPS_v1','fixture_sha256':digest(data),'state_sha256':digest(cert),
         'endpoint_Ha':str(u),'norm_enclosure':list(map(lambda v:str(F(v)),state.norm)),
         'endpoint_inequality_enclosure':list(map(lambda v:str(F(v)),gap)),
         'transfer_steps':state.steps,'enclosure_seconds':time.monotonic()-start,
         'arithmetic':'IEEE-754 binary64; adjacent outward bounds after every addition and multiplication; fastmath=False',
         'status':'certified_upper' if gap[0]>=0 else 'endpoint_refuted' if gap[1]<0 else 'ambiguous_refused'}
    if rec['status']=='ambiguous_refused' and fallback:
        exact=exact_check(data,cert)
        rec['exact_fallback']=exact
        rec['status']='certified_upper_exact_fallback' if F(exact['upper_Ha'])<=u else 'endpoint_refuted_exact_fallback'
    rec['total_seconds']=time.monotonic()-start
    return rec

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('endpoint');p.add_argument('output');p.add_argument('--fallback',action='store_true');a=p.parse_args()
    rec=check(json.loads(Path(a.fixture).read_text()),json.loads(Path(a.state).read_text()),a.endpoint,a.fallback)
    Path(a.output).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
