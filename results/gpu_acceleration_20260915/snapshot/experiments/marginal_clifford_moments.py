"""Exact moment contractions in a local valence-dimer Clifford frame.

Raw Pauli convention X^x Z^z keeps every Jordan-Wigner coefficient rational.
A creation operator is [X_j Z_<j + X_j Z_<=j]/2; annihilation has a minus.
Each four-mode block prepares (|1001>-|0110>)/sqrt(2) from frame vacuum.
This is a change of basis, with explicit finite work caps, not a scaling proof.
"""
from fractions import Fraction as F
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_symmetry_moments import simple


def jordan_wigner(h):
    result={}
    for word,coefficient in h.items():
        terms={(0,0):coefficient}
        for creation,mode in word:
            bit=1<<mode;prefix=bit-1;next_terms={}
            for (x,z),a in terms.items():
                for zz,sign in [(prefix,1),(prefix|bit,1 if creation else -1)]:
                    key=(x^bit,z^zz);value=a*F(sign*(-1)**((z&bit).bit_count()),2)
                    next_terms[key]=next_terms.get(key,0)+value
            terms={key:a for key,a in next_terms.items() if a}
        for key,a in terms.items():result[key]=result.get(key,0)+a
    result={key:simple(a) for key,a in result.items() if a}
    if any((x&z).bit_count()%2 for x,z in result):
        raise ValueError('Real Hermitian Hamiltonian must have even raw-Pauli Y parity')
    return result


def conjugate(term,gate):
    x,z=term;kind=gate[0];q=gate[1];bit=1<<q;sign=1
    if kind=='X':sign=(-1)**bool(z&bit)
    elif kind=='Z':sign=(-1)**bool(x&bit)
    elif kind=='H':
        a,b=bool(x&bit),bool(z&bit);sign=(-1)**(a and b)
        if a!=b:x^=bit;z^=bit
    elif kind=='CX':
        target=1<<gate[2]
        if x&bit:x^=target
        if z&target:z^=bit
    else:raise ValueError('Unsupported Clifford gate')
    return (x,z),sign


class CliffordMomentOracle:
    def __init__(self,certificate):
        self.base=DeterminantOracle(certificate);self.modes=self.base.modes;self.sites=self.modes//2
        if self.modes%4 or self.base.particles!=self.sites:
            raise ValueError('Dimer frame requires half filling and an even number of sites')
        if len(self.base.h)>4096:raise ValueError('Hamiltonian word budget exhausted')
        self.gates=[]
        for q in range(0,self.modes,4):
            self.gates += [('H',q),('CX',q,q+1),('CX',q,q+2),('CX',q,q+3),('Z',q),('X',q),('X',q+3)]
        terms=jordan_wigner(self.base.h);self.original_pauli_count=len(terms);self.h={}
        for term,a in terms.items():
            for gate in reversed(self.gates):term,sign=conjugate(term,gate);a*=sign
            self.h[term]=self.h.get(term,0)+a
        self.h={p:simple(a) for p,a in self.h.items() if a}
        if len(self.h)>65536:raise ValueError('Pauli term budget exhausted')
        if any((x&z).bit_count()%2 for x,z in self.h):raise ValueError('Frame Hamiltonian lost real Hermiticity')
        self.cache={};self.referenced=set()

    def _check_vector(self,v,physical=False):
        if type(v) is not dict or not 1<=len(v)<=65536:raise ValueError('Bounded nonempty sparse vector required')
        for s,a in v.items():
            if type(s) is not int or not 0<=s<1<<self.modes or (physical and not self.base.valid_state(s)):
                raise ValueError('Invalid vector coordinate or particle number')
            if (type(a) is not int and not isinstance(a,F)) or abs(a)>10**30:
                raise ValueError('Bounded exact vector coefficients required')

    def _change_basis(self,v,inverse):
        if self.sites%4:
            raise ValueError('Direct rational vector conversion requires an even number of Hadamards; use normalized_dimer_moments for odd dimer count')
        self._check_vector(v);v=dict(v)
        for gate in (list(reversed(self.gates)) if inverse else self.gates):
            kind,q=gate[:2];bit=1<<q;out={}
            for s,a in v.items():
                if kind=='H':
                    zero=s&~bit
                    out[zero]=out.get(zero,0)+a
                    out[zero|bit]=out.get(zero|bit,0)+(-a if s&bit else a)
                else:
                    t=s;sign=1
                    if kind=='X':t^=bit
                    elif kind=='Z':sign=-1 if s&bit else 1
                    elif kind=='CX' and s&bit:t^=1<<gate[2]
                    out[t]=out.get(t,0)+sign*a
            if len(out)>65536:raise ValueError('Clifford transformation support budget exhausted')
            v={s:simple(a) for s,a in out.items() if a}
        normalization=1<<(self.sites//4)
        return {s:simple(F(a,normalization)) for s,a in v.items()}

    def to_frame(self,v):
        self._check_vector(v,physical=True)
        return self._change_basis(v,True)

    def action(self,s):
        if type(s) is not int or not 0<=s<1<<self.modes:raise ValueError('Invalid frame coordinate')
        if s not in self.cache:
            if len(self.cache)>=4096:raise ValueError('Clifford frame action budget exhausted')
            result={}
            for (x,z),a in self.h.items():
                t=s^x;result[t]=result.get(t,0)+(-a if (s&z).bit_count()%2 else a)
            self.cache[s]={t:simple(a) for t,a in result.items() if a};self.referenced.update(self.cache[s]);self.referenced.add(s)
        return self.cache[s]

    def apply(self,v):
        out={}
        for s,a in v.items():
            for t,b in self.action(s).items():
                if t not in out and len(out)>=65536:raise ValueError('Frame vector support budget exhausted')
                out[t]=out.get(t,0)+a*b
        return {s:simple(a) for s,a in out.items() if a}

    def moments(self,vectors,maximum_order):
        if type(vectors) is not list or not 1<=len(vectors)<=32:raise ValueError('One to32 physical boundaries required')
        return self._moments([self.to_frame(v) for v in vectors],maximum_order)

    def dimer_moments(self,maximum_order):
        # The normalized frame vacuum maps to the normalized dimer product.
        # Carry its physical norm outside the vector to avoid sqrt(2) when
        # there is an odd number of dimers.
        M,receipt=self.normalized_dimer_moments(maximum_order)
        norm=1<<(self.sites//2)
        return [[[x*norm for x in row] for row in matrix] for matrix in M],dict(receipt,boundary_norm=norm)

    def normalized_dimer_moments(self,maximum_order):
        M,receipt=self._moments([{0:1}],maximum_order)
        return M,dict(receipt,boundary_norm=1,normalized_dimer_product=True)

    def _moments(self,boundaries,maximum_order):
        if type(maximum_order) is not int or not 0<=maximum_order<=24:raise ValueError('Moment order must lie in0..24')
        powers=[boundaries]
        for _ in range((maximum_order+1)//2):powers.append([self.apply(v) for v in powers[-1]])
        moments=[]
        for k in range(maximum_order+1):
            left,right=powers[k//2],powers[k-k//2]
            m=[[sum(a*v.get(s,0) for s,a in u.items()) for v in right] for u in left]
            if any(m[i][j]!=m[j][i] for i in range(len(m)) for j in range(len(m))):raise ValueError('Exact frame moment Hermiticity failed')
            moments.append(m)
        return moments,{'maximum_moment_order':maximum_order,'maximum_power':len(powers)-1,'boundary_count':len(boundaries),
            'initial_frame_supports':list(map(len,boundaries)),'maximum_frame_support':max(len(v) for layer in powers for v in layer),
            'source_frame_states':len(self.cache),'referenced_frame_states':len(self.referenced),
            'original_pauli_terms':self.original_pauli_count,'frame_pauli_terms':len(self.h),
            'scope':'Exact local Clifford change of basis and sparse frame moments. Frame coordinates are not fixed-N determinants; the supplied physical boundaries and number-conserving H preserve the physical state. Fixed4096 frame-source and65536 vector-support caps. No scalable closure claim.'}
