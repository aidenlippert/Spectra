"""Exact charge-flow MPS norm and CAR expectations, without a determinant basis.

Tensor entries are integers divided by one common denominator. Its power
cancels from the Rayleigh quotient; every transfer uses Python integers.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import lcm
import argparse
import json
from pathlib import Path
import sys
import time

from experiments.marginal_symbolic import decode, hermitian, validate_word
from research.molecular_collective_20260913.core import digest

I=(1,0,0,1); Z=(1,0,0,-1); CP=(0,0,1,0); AN=(0,1,0,0)


def mm(a,b):
    return (a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3])


def local_word(word,modes):
    factors=[I]*modes
    for c,j in validate_word(word,modes,2*modes):
        for k in range(j): factors[k]=mm(factors[k],Z)
        factors[j]=mm(factors[j],CP if c else AN)
    return tuple(factors)


class State:
    def __init__(self,data,cert):
        keys={'kind','fixture_sha256','modes','particles','spin_counts','denominator','bond_charges','tensors'}
        if set(cert)!=keys or cert['kind']!='integer_charge_mps_v1': raise ValueError('MPS schema')
        if cert['fixture_sha256']!=digest(data): raise ValueError('Fixture binding')
        m=cert['modes']; n=cert['particles']; den=cert['denominator']
        if type(m) is not int or m<1 or (m,n)!=(data['modes'],data['particles']): raise ValueError('Sector binding')
        if type(n) is not int or type(den) is not int or den<=0: raise ValueError('Integer sector and denominator')
        target=cert['spin_counts']
        if not isinstance(target,list) or len(target)!=2 or any(type(x) is not int or x<0 for x in target) or sum(target)!=n: raise ValueError('Spin counts')
        charges=cert['bond_charges']; tensors=cert['tensors']
        if len(charges)!=m+1 or len(tensors)!=m or charges[0]!=[[0,0]] or charges[-1]!=[target]: raise ValueError('Scalar charge boundaries')
        for cut, qs in enumerate(charges):
            if not qs: raise ValueError('Empty bond')
            for q in qs:
                if not isinstance(q,list) or len(q)!=2 or any(type(v) is not int or v<0 for v in q): raise ValueError('Integer charges')
                if q[0]>min(target[0],(cut+1)//2) or q[1]>min(target[1],cut//2): raise ValueError('Unreachable charge')
        self.rows=[]; entries=0
        for i,edges in enumerate(tensors):
            rows=[defaultdict(list),defaultdict(list)]; seen=set()
            for edge in edges:
                if len(edge)!=4 or any(type(v) is not int for v in edge): raise ValueError('Integer tensor entries')
                a,s,b,v=edge
                if s not in (0,1) or not 0<=a<len(charges[i]) or not 0<=b<len(charges[i+1]) or not v: raise ValueError('Tensor index or zero entry')
                if (a,s,b) in seen: raise ValueError('Duplicate tensor entry')
                seen.add((a,s,b)); q=list(charges[i][a]); q[i%2]+=s
                if q!=charges[i+1][b]: raise ValueError('Charge flow violation')
                rows[s][a].append((b,v)); entries+=1
            self.rows.append(rows)
        self.m=m; self.den=den; self.cert=cert
        self.stats={'bond_dimensions':list(map(len,charges)),'tensor_nonzero_entries':entries,
                    'maximum_transfer_entries':0,'transfer_steps':0,'integer_multiplications':0,
                    'enumerated_determinants':0,'many_body_matrix_entries':0}
        self.norm_integer=self.contract_factors([I]*m)
        if self.norm_integer<=0: raise ValueError('Nonpositive MPS norm')

    def step(self,E,site,O):
        rows=self.rows[site]; out=defaultdict(int); work=0
        for s in (0,1):
            if not any(O[2*s+t] for t in (0,1)): continue
            mid=defaultdict(int)
            for (a,b),ev in E.items():
                for r,av in rows[s].get(a,()): mid[r,b]+=av*ev; work+=1
            for t in (0,1):
                ov=O[2*s+t]
                if not ov: continue
                for (r,b),v in mid.items():
                    if not v: continue
                    for u,bv in rows[t].get(b,()): out[r,u]+=ov*v*bv; work+=1
        E={k:v for k,v in out.items() if v}
        self.stats['transfer_steps']+=1;self.stats['integer_multiplications']+=work
        self.stats['maximum_transfer_entries']=max(self.stats['maximum_transfer_entries'],len(E))
        return E

    def contract_factors(self,factors):
        if len(factors)!=self.m: raise ValueError('Local operator length')
        E={(0,0):1}
        for i,O in enumerate(factors): E=self.step(E,i,O)
        return E.get((0,0),0)

    def moment(self,word):
        return F(self.contract_factors(local_word(word,self.m)),self.norm_integer)

    def expectation(self,poly):
        """Shared-prefix trie; only one branch of transfer environments is live."""
        hden=lcm(*(F(c).denominator for c in poly.values())) if poly else 1
        trie={}
        for word,c in poly.items():
            node=trie
            for op in local_word(word,self.m): node=node.setdefault(op,{})
            node[None]=node.get(None,0)+int(c*hden)
        def walk(node,E,site):
            if site==self.m: return node.get(None,0)*E.get((0,0),0)
            return sum(walk(child,self.step(E,site,op),site+1) for op,child in node.items() if op is not None)
        return F(walk(trie,{(0,0):1},0),hden*self.norm_integer)


def check(data,cert):
    start=time.monotonic();state=State(data,cert)
    h=decode(data['hamiltonian'],data['modes'],4)
    if not hermitian(h) or any(sum(2*c-1 for c,_ in w) for w in h): raise ValueError('Hermitian number-conserving H required')
    upper=state.expectation(h)
    return {'status':'certified_upper','fixture_sha256':digest(data),'state_sha256':digest(cert),
            'upper_Ha':str(upper),'upper_float_Ha':float(upper),
            'norm_squared':str(F(state.norm_integer,state.den**(2*state.m))),
            'particles':data['particles'],'spin_counts':cert['spin_counts'],'total_spin_certified':False,
            'denominator':state.den,'stats':state.stats,'replay_seconds':time.monotonic()-start}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('--output');a=p.parse_args()
    result=check(json.loads(Path(a.fixture).read_text()),json.loads(Path(a.state).read_text()))
    if any(x in sys.modules for x in ('numpy','scipy','quimb','pyscf')): raise AssertionError('Numerical import in exact accepting path')
    if a.output: Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
