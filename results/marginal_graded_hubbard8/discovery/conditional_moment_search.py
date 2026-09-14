"""Physical conditional-moment coordinates for fixed-metric LP discovery.

Counts use binomial formulas; no configurations are enumerated. Any loss of
information or numerical conditioning is caught by full polynomial export.
"""
from pathlib import Path
from math import comb
from unittest.mock import patch
import json,time
import numpy as np
from scipy.sparse import csc_matrix
from experiments.marginal_moment_pricing import MomentDictionary
from experiments.marginal_joint_spinflip import flip_mask,flip_label,_orbits
from experiments.marginal_fixed_metric_pricing import solve_part

root=Path('results/marginal_graded_hubbard8');candidate=json.loads((root/'quadratic_metric_candidate.json').read_text());out=root/'fixed_conditional'
original_init=MomentDictionary.__init__;original_price=MomentDictionary.price

def initialized(self,*args,**kwargs):
    original_init(self,*args,**kwargs)
    start=time.monotonic();m=self.sites;p=self.ring.target;even=np.uint64(self.ring.spin_masks[0]);odd=np.uint64(self.ring.spin_masks[1])
    choose=np.zeros((m+1,m+1),dtype=np.int64)
    for n in range(m+1):
        for k in range(n+1):choose[n,k]=comb(n,k)
    def C(n,k):
        valid=(n>=0)&(k>=0)&(k<=n)
        return np.where(valid,choose[np.clip(n,0,m),np.clip(k,0,m)],0)
    def event(required,occupied):
        required=np.asarray(required,dtype=np.uint64);occupied=np.asarray(occupied,dtype=np.uint64)
        na=np.bitwise_count(required&even).astype(int);nb=np.bitwise_count(required&odd).astype(int)
        a=np.bitwise_count(occupied&even).astype(int);b=np.bitwise_count(occupied&odd).astype(int)
        full=C(m-na,p-a)*C(m-nb,p-b)
        empty=required^occupied
        up=(occupied&even)|((empty&odd)>>np.uint64(1));down=((occupied&odd)>>np.uint64(1))|(empty&even)
        forced=np.bitwise_count(up|down).astype(int);au=np.bitwise_count(up).astype(int)
        valence=np.where((up&down)==0,C(m-forced,p-au),0)
        return full-valence
    basis=np.array(self.quotient.basis,dtype=np.uint64);qp=[self.quotient.index[flip_mask(int(mask),m)] for mask in basis];orbits=_orbits(qp);indices=np.zeros(len(basis),dtype=int)
    for j,orbit in enumerate(orbits):indices[orbit]=j
    representatives=np.array([basis[orbit[0]] for orbit in orbits],dtype=np.uint64);counts=event(representatives,representatives)
    if np.any(counts<=0):raise ValueError('Conditional test monomial has no Q completion')
    T=np.empty((self.qrows,self.qrows))
    for row,support in enumerate(representatives):
        union=support|basis
        T[row]=np.bincount(indices,weights=event(union,union),minlength=self.qrows)/counts[row]
    # Check sampled atom columns against direct binomial event counts.
    maximum=0.;labels=self.labels();sample=labels[::max(1,len(labels)//40)]
    for label in sample:
        orbit=sorted({tuple(label),tuple(flip_label(label,m))})
        expected=np.zeros(self.qrows)
        for family,required,occupied in orbit:
            R=np.uint64(required);O=np.uint64(occupied);empty=R^O;conflict=(representatives&empty)!=0
            rr=representatives|R;oo=representatives|O
            value=np.where(conflict,0,event(rr,oo)).astype(float)
            if family=='charge':
                value=-value
                for i in range(m):
                    pair=np.uint64(3 << (2*i))
                    value+=np.where(conflict|((pair&empty)!=0),0,event(rr|pair,oo|pair))
            expected+=value/(counts*len(orbit))
        actual=T@self.column(label)
        maximum=max(maximum,float(max(abs(actual-expected))))
    if maximum>1e-8:raise ValueError('Conditional-moment coordinate check failed')
    self.conditional_T=T;self.projection=csc_matrix(T)@self.projection;self.one=T@self.one
    self.stats.update(conditional_rows=self.qrows,conditional_transform_entries=T.size,
        conditional_construction_seconds=time.monotonic()-start,checked_columns=len(sample),maximum_direct_count_error=maximum,
        conditional_scope='Exact binomial physical-Q test functions converted to floating proposal coordinates; original polynomial export is unchanged. Invertibility is not used as an acceptance shortcut.')

def price(self,dual,active,batch):
    transformed=np.array(dual,copy=True);r=self.qrows
    transformed[:r]=self.conditional_T.T@transformed[:r];transformed[r:2*r]=self.conditional_T.T@transformed[r:2*r]
    return original_price(self,transformed,active,batch)

with patch.object(MomentDictionary,'__init__',initialized),patch.object(MomentDictionary,'price',price),patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No states')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full lift')):
    result=solve_part(candidate,'numerator',out/'numerator',native=True,time_limit=240)
print(json.dumps(result),flush=True)
