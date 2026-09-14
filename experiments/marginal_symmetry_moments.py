"""Exact full moments in signed spin/particle-hole orbits.

The orbit basis is unnormalized: E_r=sum_t phase_r(t)|t>. Its norm squared
is orbit size. Only one determinant CAR action per surviving orbit is used.
This is a finite-group compression, not a polynomial many-body scaling claim.
"""
from fractions import Fraction as F
from math import comb
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_symbolic import canonical, add, scale, number_shift
from experiments.marginal_enlarged_schur import solve_positive


def simple(x):
    return x.numerator if isinstance(x,F) and x.denominator==1 else x


def projected_moments(moments):
    """Recover V^T H Q(QHQ)^k Q H V from full nonorthogonal moments.

    These inputs must come from a verified contraction. Matrix checks alone
    do not establish that arbitrary submitted moments belong to a physical H.
    """
    if not 3<=len(moments)<=25:raise ValueError('Three to25 full moments required')
    n=len(moments[0])
    if not 1<=n<=32:raise ValueError('Moment matrix dimension exceeds32')
    for m in moments:
        if len(m)!=n or any(len(row)!=n for row in m):raise ValueError('Inconsistent moment dimensions')
        if any(type(x) is not int and not isinstance(x,F) for row in m for x in row):
            raise ValueError('Exact integer or rational moments required')
        if any(m[i][j]!=m[j][i] for i in range(n) for j in range(n)):raise ValueError('Hermitian moments required')
    def mul(a,b):return [[sum(x*b[k][j] for k,x in enumerate(row)) for j in range(n)] for row in a]
    G=moments[0];A=[None]+[solve_positive(G,m) for m in moments[1:]]
    C=[[[F(i==j) for j in range(n)] for i in range(n)]]
    for k in range(1,len(moments)):
        products=[mul(A[j],C[k-j]) for j in range(1,k+1)]
        C.append([[-sum(p[i][j] for p in products) for j in range(n)] for i in range(n)])
    result=[[[-simple(x) for x in row] for row in mul(G,c)] for c in C[2:]]
    if any(m[i][j]!=m[j][i] for m in result for i in range(n) for j in range(n)):
        raise ValueError('Projected moment Hermiticity failed')
    return result


class SymmetryMomentOracle:
    def __init__(self, certificate):
        self.base=DeterminantOracle(certificate)
        self.modes=self.base.modes;self.sites=self.modes//2
        if self.modes%8 or self.base.particles!=self.sites:
            raise ValueError('Orbit convention requires half filling and a multiple of four sites')
        self.alpha_mask=sum(1<<(2*i) for i in range(self.sites))
        h=self.base.h
        if any(sum((2*c-1) for c,i in w if i%2==0) for w in h):
            raise ValueError('Spin population conservation required')
        flip={tuple((c,i^1) for c,i in w):a for w,a in h.items()}
        if canonical(flip)!=h:
            raise ValueError('Exact spin-exchange Hamiltonian symmetry required')
        ph={tuple((1-c,i) for c,i in w):a*(-1)**sum(i//2 for _,i in w) for w,a in h.items()}
        difference=add(canonical(ph),scale(h,-1))
        multiplier=-difference.get((),F(0))/self.base.particles
        if difference!=scale(number_shift(self.modes,self.base.particles),multiplier):
            raise ValueError('Particle-hole symmetry modulo fixed number required')
        self.symmetry_receipt={'spin_exchange_exact':True,'particle_hole_number_multiplier':str(multiplier),
                               'particle_hole_fixed_number_exact':True}
        self.orbits={};self.cache={}

    def valid_state(self,s):
        return self.base.valid_state(s) and (s&self.alpha_mask).bit_count()==self.sites//2

    def flip(self,s):
        swapped=((s&self.alpha_mask)<<1)|((s>>1)&self.alpha_mask)
        doublons=sum(((s>>(2*i))&3)==3 for i in range(self.sites))
        return swapped,(-1)**doublons

    def particle_hole(self,s):
        phase=(-1)**sum(i//2+i for i in range(self.modes) if s&(1<<i))
        return self.base.all_bits^s,phase

    def orbit(self,s):
        if not self.valid_state(s):raise ValueError('Balanced physical determinant required')
        if s in self.orbits:return self.orbits[s]
        f,sf=self.flip(s);c,sc=self.particle_hole(s);fc,sfc=self.flip(c)
        raw={};zero=False
        for t,phase in [(s,1),(f,sf),(c,sc),(fc,sc*sfc)]:
            if t in raw and raw[t]!=phase:zero=True
            raw[t]=phase
        if zero:
            for t in raw:self.orbits[t]=None
            return None
        r=min(raw);sign=raw[r];phases={t:phase*sign for t,phase in raw.items()}
        for t in raw:self.orbits[t]=(r,phases[t],len(raw),phases)
        return self.orbits[s]

    def compress(self,v):
        if type(v) is not dict or not 1<=len(v)<=4096:
            raise ValueError('One to4096 sparse boundary coefficients required')
        if any((type(a) is not int and not isinstance(a,F)) or abs(a)>10**30 for a in v.values()):
            raise ValueError('Bounded exact integer or rational boundary coefficients required')
        result={};checked=set()
        for s,a in v.items():
            if not a:continue
            entry=self.orbit(s)
            if entry is None:raise ValueError('Nonzero vector on a forbidden signed orbit')
            r,phase,size,orbit=entry
            if r in checked:continue
            coefficient=v.get(r,0)
            if any(v.get(t,0)!=coefficient*p for t,p in orbit.items()):
                raise ValueError('Input vector is not invariant under both signed symmetries')
            checked.add(r);result[r]=simple(coefficient)
        return result

    def action(self,r):
        entry=self.orbit(r)
        if entry is None or entry[0]!=r:raise ValueError('Canonical surviving orbit representative required')
        if r not in self.cache:
            size=entry[2];result={}
            for t,a in self.base.action(r).items():
                target=self.orbit(t)
                if target is None:continue
                u,phase,target_size,_=target
                result[u]=result.get(u,0)+F(size,target_size)*phase*a
            self.cache[r]={u:simple(a) for u,a in result.items() if a}
        return self.cache[r]

    def apply(self,v):
        out={}
        for r,a in v.items():
            for t,b in self.action(r).items():out[t]=out.get(t,0)+a*b
        return {r:simple(a) for r,a in out.items() if a}

    def dot(self,a,b):
        return sum(x*b.get(r,0)*self.orbit(r)[2] for r,x in a.items())

    def moments(self,vectors,maximum_order):
        if type(maximum_order) is not int or not 0<=maximum_order<=24:
            raise ValueError('Moment order must lie in0..24')
        if type(vectors) is not list or not 1<=len(vectors)<=32:raise ValueError('One to32 boundary vectors required')
        powers=[[self.compress(v) for v in vectors]]
        for _ in range((maximum_order+1)//2):powers.append([self.apply(v) for v in powers[-1]])
        moments=[]
        for k in range(maximum_order+1):
            left,right=powers[k//2],powers[k-k//2]
            m=[[self.dot(a,b) for b in right] for a in left]
            if any(m[i][j]!=m[j][i] for i in range(len(m)) for j in range(len(m))):
                raise ValueError('Exact moment Hermiticity failed')
            moments.append(m)
        return moments,{'maximum_moment_order':maximum_order,'maximum_power':len(powers)-1,
            'boundary_count':len(vectors),'maximum_orbit_support':max(len(v) for layer in powers for v in layer),
            'unique_determinant_source_actions':len(self.base.cache),'referenced_determinants':self.base.referenced_state_count(),
            'orbit_cache_states':len(self.orbits),'balanced_sector_dimension':comb(self.sites,self.sites//2)**2,
            'materialized_full_sector':False,'symmetry_checks':self.symmetry_receipt,
            'scope':'Exact signed finite-group orbit contraction from prescribed invariant boundary vectors. Orbits and their images are generated on demand; no full sector matrix. Orbit count can still grow exponentially.'}
