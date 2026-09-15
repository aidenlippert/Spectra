"""Numerical MPS moment proposals, never an accepting path or a state vector."""
import json
import time
from pathlib import Path
import numpy as np
from numba import njit
from experiments.marginal_symbolic import word_product
from experiments.marginal_coefficient import dagger
from research.correlated_pair_20260913.mps_exact import State, I, local_word

@njit(cache=True)
def transfer(E,edges,O,right):
    left=E.shape[0]
    result=np.zeros((right,right))
    for s in range(2):
        if O[2*s]==0 and O[2*s+1]==0: continue
        mid=np.zeros((right,left))
        for k in range(len(edges)):
            a,p,b,v=edges[k]
            if int(p)!=s: continue
            for j in range(left):
                mid[int(b),j]+=v*E[int(a),j]
        for k in range(len(edges)):
            a,p,b,v=edges[k]; ov=O[2*s+int(p)]
            if ov==0: continue
            for j in range(right):
                result[j,int(b)]+=ov*mid[j,int(a)]*v
    return result

class Proposal(State):
    def contract_factors(self,factors):
        # State.__init__ performs every existing structural gate. This sentinel
        # only skips its expensive *exact* norm in this proposal-only class.
        # A numerical norm is computed and checked below, never certified here.
        return 1

    def __init__(self,data,cert):
        super().__init__(data,cert)
        self.edges=[np.array([[a,s,b,v/self.den] for a,s,b,v in site],dtype=float).reshape(-1,4) for site in cert['tensors']]
        self.dims=list(map(len,cert['bond_charges']))
        self.cache={};self.steps=0
        self.norm=self.evaluate_factors([I]*self.m)
        if not np.isfinite(self.norm) or self.norm<=0: raise ValueError('Numerical proposal norm')

    def evaluate_factors(self,factors):
        E=np.ones((1,1))
        for site,op in enumerate(factors):
            E=transfer(E,self.edges[site],np.array(op,dtype=np.int64),self.dims[site+1]);self.steps+=1
        return float(E[0,0])

    def fill(self,words):
        trie={};todo=[w for w in words if w not in self.cache]
        for word in todo:
            # A fixed spin-count MPS makes charge-changing moments zero.
            if any(sum(2*c-1 for c,i in word if i%2==p) for p in (0,1)):
                self.cache[word]=0.;continue
            ops=local_word(word,self.m)
            if any(not any(op) for op in ops):self.cache[word]=0.;continue
            node=trie
            for op in ops:node=node.setdefault(op,{})
            node.setdefault(None,[]).append(word)
        def walk(node,E,site):
            if site==self.m:
                for word in node[None]:self.cache[word]=float(E[0,0])/self.norm
                return
            for op,child in node.items():
                self.steps+=1
                walk(child,transfer(E,self.edges[site],np.array(op,dtype=np.int64),self.dims[site+1]),site+1)
        if trie:walk(trie,np.ones((1,1)),0)

    def matrix(self,words):
        expansions={(i,j):word_product(dagger(left),right) for i,left in enumerate(words) for j,right in enumerate(words[:i+1])}
        self.fill({w for terms in expansions.values() for w,c in terms})
        C=np.zeros((len(words),len(words)))
        for (i,j),terms in expansions.items():
            C[i,j]=C[j,i]=sum(float(c)*self.cache[w] for w,c in terms)
        return C

