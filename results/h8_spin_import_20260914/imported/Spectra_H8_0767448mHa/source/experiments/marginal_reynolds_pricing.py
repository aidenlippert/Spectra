"""Invariant coefficient coordinates for finite-group moment pricing.

Numerical row selection is a proposal reduction only. Export must expand every
atom orbit and pass the unchanged exact polynomial verifier.
"""
from fractions import Fraction as F
import time
from experiments.marginal_moment_pricing import MomentDictionary
from experiments.marginal_joint_spinflip import flip_label,flip_mask,_orbits
from experiments.marginal_joint_coefficient_constructor import vanishing_feature
from experiments.marginal_polynomial_metric import JointPolynomial


def reflect_mask(mask,sites):
    return sum(((mask>>(2*i))&3)<<(2*(sites-1-i)) for i in range(sites))


def complement_polynomial(polynomial):
    result={}
    for mask,coefficient in polynomial.items():
        subset=mask
        while True:
            result[subset]=result.get(subset,0)+coefficient*(-1 if subset.bit_count()%2 else 1)
            if not subset:break
            subset=(subset-1)&mask
    return {m:v for m,v in result.items() if v}


class ReynoldsMomentDictionary(MomentDictionary):
    def __init__(self,data,gamma,feature_degree=2,proof_degree=6,reflection=True,particle_hole=True,coordinates="coefficient"):
        import numpy as np
        from scipy.sparse import csc_matrix,block_diag
        from scipy.linalg import qr
        if type(reflection) is not bool or type(particle_hole) is not bool:
            raise ValueError('Boolean symmetry flags required')
        if coordinates not in ('coefficient','conditional'):raise ValueError('Supported Reynolds coordinates required')
        self.reflection,self.particle_hole=reflection,particle_hole
        super().__init__(data,gamma,feature_degree=feature_degree,proof_degree=proof_degree)
        started=time.monotonic();quotient=self.quotient
        # Each charge monomial has definite particle-hole parity. Restrict to
        # even features, then check the actual compiled target exactly.
        selected=[j for j,orbit in enumerate(self.orbits) if not particle_hole or all(sum(p)%2==0 for p in orbit)]
        self.retained_feature_indices=selected
        self.orbits=[self.orbits[j] for j in selected]
        self.metric=self.metric[:,selected]
        for orbit in self.orbits:
            ring=JointPolynomial(dict(data,polynomial_metric=vanishing_feature(orbit,self.sites)))
            weight,numerator,_=ring.compile(F(gamma))
            for polynomial in (weight,numerator):self.check_invariant(polynomial)
        qp=[quotient.index[flip_mask(m,self.sites)] for m in quotient.basis]
        spin_orbits=_orbits(qp);ri=[];ci=[];values=[]
        group_size=(2 if reflection else 1)*(2 if particle_hole else 1)
        for j,members in enumerate(spin_orbits):
            image={}
            for member in members:
                mask=quotient.basis[member]
                supports=(mask,reflect_mask(mask,self.sites)) if reflection else (mask,)
                for support in supports:
                    image[support]=image.get(support,0)+1
                    if particle_hole:
                        for subset,value in complement_polynomial({support:1}).items():image[subset]=image.get(subset,0)+value
            for mask,value in image.items():
                if value:ri.append(self.index[mask]);ci.append(j);values.append(value/group_size)
        images=csc_matrix((values,(ri,ci)),shape=(len(self.monomials),self.qrows))
        S=(self.moment_map@images).toarray();rank=round(float(np.trace(S)))
        if np.max(abs(S@S-S))>1e-8:raise ValueError('Reynolds idempotence diagnostic failed')
        row_map=S
        if coordinates=='conditional':
            from experiments.marginal_conditional_coordinates import conditional_map
            T,conditional_stats=conditional_map(self)
            row_map=T@S;self.stats.update(conditional_stats)
        _,R,pivots=qr(row_map.T,mode='economic',pivoting=True,check_finite=False)
        if not 0<rank<=self.qrows or abs(R[rank-1,rank-1])<1e-8 or (rank<self.qrows and np.max(abs(np.diag(R)[rank:]))>1e-8):
            raise ValueError('Reynolds rank diagnostic failed')
        self.original_rows=self.qrows;self.reynolds_L=row_map[pivots[:rank]];L=csc_matrix(self.reynolds_L)
        self.projection=L@self.projection;self.one=self.reynolds_L@self.one
        self.metric=block_diag((L,L,csc_matrix([[1.]])),format='csc')@self.metric;self.qrows=rank
        self.stats.update(reynolds_rows=rank,reynolds_coordinates=coordinates,reynolds_original_rows=self.original_rows,reynolds_metric_features=len(self.orbits),reynolds_seconds=time.monotonic()-started,
            reynolds_actions=[name for enabled,name in ((reflection,'reflection'),(particle_hole,'particle-hole')) if enabled],
            reynolds_scope='Exact quotient invariance of retained compiled metric and numerator columns; numerical independent-row proposals; group expansion and full original exact replay required.')

    def check_invariant(self,polynomial):
        target=self.quotient.normal(polynomial)
        if self.reflection and self.quotient.normal({reflect_mask(m,self.sites):v for m,v in polynomial.items()})!=target:
            raise ValueError('Reflection target invariance failed')
        if self.particle_hole and self.quotient.normal(complement_polynomial(polynomial))!=target:
            raise ValueError('Particle-hole target invariance failed')

    def atom_members(self,label):
        family,R,O=label
        if O & ~R:raise ValueError('Occupied mask must lie inside required support')
        items={(family,R,O)}
        if self.particle_hole:items|={(f,r,r^o) for f,r,o in tuple(items)}
        if self.reflection:items|={(f,reflect_mask(r,self.sites),reflect_mask(o,self.sites)) for f,r,o in tuple(items)}
        return tuple(sorted({min(member,tuple(flip_label(member,self.sites))) for member in items}))

    def labels(self,degree=6):
        return sorted({min(self.atom_members(label)) for label in super().labels(degree)})

    def price(self,dual,active,batch):
        import numpy as np
        r,old=self.qrows,self.original_rows;L=self.reynolds_L
        lifted=np.r_[L.T@dual[:r],L.T@dual[r:2*r],0.]
        expanded={(b,member) for b,label in active for member in self.atom_members(label)}
        self.qrows=old
        try:raw,maximum,checked=super().price(lifted,expanded,4*batch)
        finally:self.qrows=r
        chosen=[];seen=set(active);counts=[0,0]
        for b,label in raw:
            canonical=min(self.atom_members(label))
            if (b,canonical) not in seen and counts[b]<batch:
                chosen.append((b,canonical));seen.add((b,canonical));counts[b]+=1
        return chosen,maximum,checked
