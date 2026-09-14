"""Reflection/particle-hole averaging in a fixed-metric coefficient search.

Exact target invariance is checked; numerical row selection only proposes a
certificate. Expanded atom weights must pass the full original exact replay.
"""
from pathlib import Path
from unittest.mock import patch
import json,time,shutil
import numpy as np
from scipy.sparse import csc_matrix
from scipy.linalg import qr
from experiments.marginal_moment_pricing import MomentDictionary
from experiments.marginal_joint_spinflip import flip_mask,flip_label,_orbits
from experiments.marginal_polynomial_metric import JointPolynomial
import experiments.marginal_fixed_metric_pricing as fixed

root=Path('results/marginal_graded_hubbard8');candidate=json.loads((root/'quadratic_metric_candidate.json').read_text());out=root/'fixed_reynolds'
original_init=MomentDictionary.__init__;original_price=MomentDictionary.price;original_labels=MomentDictionary.labels;original_export=fixed.export_part

def reflect(mask,sites=8):
    return sum(((mask>>(2*i))&3)<<(2*(sites-1-i)) for i in range(sites))

def orbit(label):
    family,R,O=label;items=set()
    for r,o in ((R,O),(R,R^O),(reflect(R),reflect(O)),(reflect(R),reflect(R^O))):
        member=(family,r,o);items.add(min(member,tuple(flip_label(member,8))))
    return sorted(items)

def complement(polynomial):
    result={}
    for mask,coefficient in polynomial.items():
        subset=mask
        while True:
            result[subset]=result.get(subset,0)+coefficient*(-1 if subset.bit_count()%2 else 1)
            if not subset:break
            subset=(subset-1)&mask
    return {m:v for m,v in result.items() if v}

def initialized(self,*args,**kwargs):
    original_init(self,*args,**kwargs)
    if self.sites!=8:raise ValueError('Eight-site prototype only')
    started=time.monotonic();quotient=self.quotient
    ring=JointPolynomial(candidate);w,k,cost=ring.compile(candidate['target_lower'])
    for polynomial in (w,k):
        target=quotient.normal(polynomial)
        if quotient.normal({reflect(m):v for m,v in polynomial.items()})!=target:raise ValueError('Reflection target invariance failed')
        if quotient.normal(complement(polynomial))!=target:raise ValueError('Particle-hole target invariance failed')
    qp=[quotient.index[flip_mask(m,8)] for m in quotient.basis];spin_orbits=_orbits(qp)
    ri=[];ci=[];values=[]
    for j,members in enumerate(spin_orbits):
        image={}
        for member in members:
            mask=quotient.basis[member]
            for support in (mask,reflect(mask)):
                image[support]=image.get(support,0)+1
                subset=support
                while True:
                    image[subset]=image.get(subset,0)+(-1 if subset.bit_count()%2 else 1)
                    if not subset:break
                    subset=(subset-1)&support
        for mask,value in image.items():
            if value:ri.append(self.index[mask]);ci.append(j);values.append(value/4)
    images=csc_matrix((values,(ri,ci)),shape=(len(self.monomials),self.qrows))
    S=(self.moment_map@images).toarray();rank=round(float(np.trace(S)))
    if np.max(abs(S@S-S))>1e-8:raise ValueError('Reynolds idempotence check failed')
    _,R,pivots=qr(S.T,mode='economic',pivoting=True,check_finite=False)
    if not 0<rank<self.qrows or abs(R[rank-1,rank-1])<1e-8 or np.max(abs(np.diag(R)[rank:]))>1e-8:raise ValueError('Reynolds rank check failed')
    L=S[pivots[:rank]];self.original_rows=self.qrows;self.reynolds_L=L
    self.projection=csc_matrix(L)@self.projection;self.one=L@self.one;self.qrows=rank
    self.stats.update(reynolds_rows=rank,reynolds_original_rows=self.original_rows,reynolds_trace=float(np.trace(S)),reynolds_seconds=time.monotonic()-started,
        reynolds_scope='Exact coefficient target invariance under reflection and particle-hole; numerical independent-row proposal reduction; all group weights expanded for exact original replay.')
    print(json.dumps({'phase':'reynolds_ready','rows':rank,'original_rows':self.original_rows,'seconds':time.monotonic()-started}),flush=True)

def labels(self,degree=6):
    return sorted({min(orbit(label)) for label in original_labels(self,degree)})

def price(self,dual,active,batch):
    r=self.qrows;old=self.original_rows;L=self.reynolds_L
    lifted=np.r_[L.T@dual[:r],L.T@dual[r:2*r],0.]
    expanded={(b,member) for b,label in active for member in orbit(label)}
    self.qrows=old
    try:raw,maximum,checked=original_price(self,lifted,expanded,4*batch)
    finally:self.qrows=r
    chosen=[];seen=set(active);counts=[0,0]
    for b,label in raw:
        canonical=min(orbit(label))
        if (b,canonical) not in seen and counts[b]<batch:
            chosen.append((b,canonical));seen.add((b,canonical));counts[b]+=1
    return chosen,maximum,checked

def export_part(ring,polynomial,scale,labels,weights,bound,denominator=10**12):
    expanded_labels=[];expanded_weights=[]
    for label,weight in zip(labels,weights):
        members=orbit(label)
        for member in members:expanded_labels.append(member);expanded_weights.append(float(weight)/len(members))
    return original_export(ring,polynomial,scale,expanded_labels,expanded_weights,bound,denominator)

with patch.object(MomentDictionary,'__init__',initialized),patch.object(MomentDictionary,'labels',labels),patch.object(MomentDictionary,'price',price),patch.object(fixed,'export_part',export_part),patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No states')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full lift')):
    result=fixed.solve_part(candidate,'numerator',out/'numerator',native=True,time_limit=240,batch=128)
if result['exact_accepted']:
    (out/'weight').mkdir(exist_ok=True)
    for filename in ('receipt.json','proof.json'):shutil.copyfile(root/'fixed_metric_native/weight'/filename,out/'weight'/filename)
    print(json.dumps(fixed.combine(candidate,out)),flush=True)
