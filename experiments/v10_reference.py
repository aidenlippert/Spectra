"""Bounded exact Duhamel construction around onsite Z + local depolarization.

This is a supplied conventional calculation, not a learned method. The proposer
exports every D_k. The checker verifies differential identities independently
of the integration routine. Frequencies and all polynomial coefficients are
complex rationals, including exactly coincident frequencies.
"""
from fractions import Fraction as F
from math import factorial
from .v7_certificate import Generator, clean, rational, BudgetExceeded

ZERO = (F(0), F(0))
ONE = (F(1), F(0))
IMAG = (F(0), F(1))


def cq(z):
    if not isinstance(z, tuple) or len(z) != 2:
        raise ValueError('complex rational pair required')
    return rational(z[0]), rational(z[1])


def add(a, b):
    return rational(a[0]+b[0]), rational(a[1]+b[1])


def neg(a):
    return -a[0], -a[1]


def mul(a, b):
    return rational(a[0]*b[0]-a[1]*b[1]), rational(a[0]*b[1]+a[1]*b[0])


def scale(a, c):
    return rational(a[0]*c), rational(a[1]*c)


def div(a, b):
    s = rational(b[0]**2+b[1]**2)
    if not s: raise ValueError('zero complex denominator')
    return rational((a[0]*b[0]+a[1]*b[1])/s), rational((a[1]*b[0]-a[0]*b[1])/s)


def conj(a):
    return a[0], -a[1]


def accumulate(out, key, value):
    z = add(out.get(key, ZERO), value)
    if z != ZERO: out[key] = z
    elif key in out: del out[key]


PAULI_TO_EIGEN = {
    'I': {'I': ONE}, 'Z': {'Z': ONE},
    'X': {'P': ONE, 'M': ONE}, 'Y': {'P': neg(IMAG), 'M': IMAG},
}
EIGEN_TO_PAULI = {
    'I': {'I': ONE}, 'Z': {'Z': ONE},
    'P': {'X': scale(ONE,F(1,2)), 'Y': scale(IMAG,F(1,2))},
    'M': {'X': scale(ONE,F(1,2)), 'Y': scale(IMAG,F(-1,2))},
}
LOCAL_PRODUCT = {
    ('Z','Z'): {'I': ONE}, ('Z','P'): {'P': ONE}, ('P','Z'): {'P': neg(ONE)},
    ('Z','M'): {'M': neg(ONE)}, ('M','Z'): {'M': ONE},
    ('P','M'): {'I': scale(ONE,F(1,2)), 'Z': scale(ONE,F(1,2))},
    ('M','P'): {'I': scale(ONE,F(1,2)), 'Z': scale(ONE,F(-1,2))},
    ('P','P'): {}, ('M','M'): {},
}


def expand_word(word, table, cap=512):
    out = {'': ONE}
    for letter in word:
        out = {p+q: mul(c,d) for p,c in out.items() for q,d in table[letter].items()}
        if len(out)>cap: raise BudgetExceeded('basis conversion support budget')
    return out


def multiply_words(left, right, cap=512):
    out = {'': ONE}
    for a,b in zip(left,right):
        local = {b:ONE} if a=='I' else ({a:ONE} if b=='I' else LOCAL_PRODUCT[a,b])
        out = {p+q:mul(c,d) for p,c in out.items() for q,d in local.items()}
        if len(out)>cap: raise BudgetExceeded('product expansion budget')
    return out


def pauli_to_eigen(op, n, cap=512):
    out = {}
    for p,c in clean(op,n,cap).items():
        for q,d in expand_word(p,PAULI_TO_EIGEN,cap).items(): accumulate(out,q,scale(d,c))
        if len(out)>cap: raise BudgetExceeded('eigenoperator support budget')
    return out


def eigen_to_pauli(op, n, cap=512):
    out = {}
    for p,c in op.items():
        if len(p)!=n or any(x not in 'IZPM' for x in p): raise ValueError('eigenoperator label')
        c = cq(c)
        for q,d in expand_word(p,EIGEN_TO_PAULI,cap).items(): accumulate(out,q,mul(c,d))
        if len(out)>cap: raise BudgetExceeded('Pauli output support budget')
    return out


class Reference:
    def __init__(self, h, gamma, n, cap=512, entries=50000):
        g = Generator(h,gamma,n,cap)
        if type(entries) is not int or not 1<=entries<=50000: raise ValueError('entry cap')
        self.h,self.gamma,self.n,self.cap,self.entries = g.h,g.gamma,n,cap,entries
        self.fields = [F(0)]*n
        self.v = {}
        for p,c in self.h.items():
            if p=='I'*n: continue  # Scalar Hamiltonian contributes no commutator.
            if p.count('Z')==1 and p.count('I')==n-1: self.fields[p.index('Z')] = c
            else: self.v[p] = c
        self.ve = pauli_to_eigen(self.v,n,cap)
        self.k_bound = rational(2*sum(map(abs,self.v.values()),F(0)))
        self.cache = {}
        self.cost = dict(column_requests=0, cache_hits=0, interaction_pairs=0,
                         complex_multiply_adds=0, integral_terms=0, repeated_frequencies=0,
                         peak_entries=0, peak_support=0)

    def frequency(self, word):
        return (-self.gamma*sum(c!='I' for c in word),
                2*sum((h if c=='P' else -h if c=='M' else F(0)) for h,c in zip(self.fields,word)))

    def column(self, word):
        self.cost['column_requests'] += 1
        if word in self.cache:
            self.cost['cache_hits'] += 1
            return self.cache[word]
        if len(self.cache)>=self.cap: raise BudgetExceeded('reference cached column budget')
        out = {}
        for p,c in self.ve.items():
            self.cost['interaction_pairs'] += 1
            for q,d in multiply_words(p,word,self.cap).items(): accumulate(out,q,mul(IMAG,mul(c,d)))
            for q,d in multiply_words(word,p,self.cap).items(): accumulate(out,q,neg(mul(IMAG,mul(c,d))))
            if len(out)>self.cap: raise BudgetExceeded('interaction column support budget')
        if len(out)>self.cap: raise BudgetExceeded('interaction column support budget')
        self.cache[word] = out
        return out

    def validate(self, terms):
        if not isinstance(terms,dict) or len(terms)>self.entries: raise BudgetExceeded('mode entry budget')
        support = set()
        for key,c in terms.items():
            if not isinstance(key,tuple) or len(key)!=4: raise ValueError('mode key')
            word,a,b,k = key
            if not isinstance(word,str) or len(word)!=self.n or any(x not in 'IZPM' for x in word): raise ValueError('mode word')
            if type(k) is not int or not 0<=k<=32: raise ValueError('mode degree')
            rational(a); rational(b); cq(c)
            if a>0: raise ValueError('positive real exponent')
            support.add(word)
        if len(support)>self.cap: raise BudgetExceeded('mode support budget')
        self.cost['peak_entries'] = max(self.cost['peak_entries'],len(terms))
        self.cost['peak_support'] = max(self.cost['peak_support'],len(support))
        return terms

    def apply_b(self, terms):
        self.validate(terms)
        out = {}
        for (p,a,b,k),c in terms.items():
            for q,d in self.column(p).items():
                accumulate(out,(q,a,b,k),mul(c,d))
                self.cost['complex_multiply_adds'] += 1
            if len(out)>self.entries: raise BudgetExceeded('interaction mode entry budget')
        return self.validate(out)

    def initial(self, op):
        return self.validate({(p,*self.frequency(p),0): c for p,c in pauli_to_eigen(op,self.n,self.cap).items()})

    def integrate(self, forcing):
        """Exact exp(mu*t) integral_0^t s^k exp((lambda-mu)*s) ds."""
        self.validate(forcing)
        out = {}
        for (p,a,b,k),c in forcing.items():
            mu = self.frequency(p)
            diff = add((a,b),neg(mu))
            if diff==ZERO:
                if k==32: raise BudgetExceeded('integrated degree budget')
                accumulate(out,(p,a,b,k+1),scale(c,F(1,k+1)))
                self.cost['repeated_frequencies'] += 1
                self.cost['integral_terms'] += 1
            else:
                # exp(lambda*t) Q(t) - exp(mu*t) Q(0), Q'+diff*Q=c*t^k.
                q = div(c,diff)
                for j in range(k,-1,-1):
                    accumulate(out,(p,a,b,j),q)
                    self.cost['integral_terms'] += 1
                    if j: q = div(scale(q,-j),diff)
                accumulate(out,(p,*mu,0),neg(q))
                self.cost['integral_terms'] += 1
            if len(out)>self.entries: raise BudgetExceeded('integrated mode entry budget')
        return self.validate(out)


def hermitian_modes(terms):
    swap = str.maketrans('PM','MP')
    for (p,a,b,k),c in terms.items():
        if terms.get((p.translate(swap),a,-b,k),ZERO)!=conj(c): return False
    return True


def derivative_minus_reference(ref, terms):
    """Checker action; does not call the proposer's integration routine."""
    out = {}
    for (p,a,b,k),c in terms.items():
        mu = ref.frequency(p)
        accumulate(out,(p,a,b,k),mul(add((a,b),neg(mu)),c))
        if k: accumulate(out,(p,a,b,k-1),scale(c,k))
    return ref.validate(out)


def at_zero(terms):
    out = {}
    for (p,a,b,k),c in terms.items():
        if k==0: accumulate(out,p,c)
    return out


def check_expansion(h, gamma, n, initial, duration, tolerance, layers, cap=512):
    """Recompute all ODE identities and bound; no trusted claimed scalar."""
    ref = Reference(h,gamma,n,cap)
    duration,tolerance = rational(duration),rational(tolerance)
    if duration<=0 or tolerance<0: raise ValueError('duration/tolerance')
    initial = clean(initial,n,cap)
    if not isinstance(layers,(tuple,list)) or not 1<=len(layers)<=25: raise ValueError('expansion order cap')
    total_entries=0
    for k,layer in enumerate(layers):
        ref.validate(layer)
        total_entries+=len(layer)
        if total_entries>50000: raise BudgetExceeded('expansion total entry budget')
        if not hermitian_modes(layer): raise ValueError('non-Hermitian expansion layer')
        if k==0:
            if layer!=ref.initial(initial): raise ValueError('reference initial layer mismatch')
        else:
            if at_zero(layer): raise ValueError('nonzero correction initial condition')
            lhs=derivative_minus_reference(ref,layer)
            rhs=ref.apply_b(layers[k-1])
            if lhs!=rhs: raise ValueError('correction differential identity mismatch')
    m = len(layers)-1
    bound=rational(sum(map(abs,initial.values()),F(0))*(ref.k_bound*duration)**(m+1)/factorial(m+1))
    return dict(status='certified' if bound<=tolerance else 'over_tolerance',
                bound=str(bound), k_bound=str(ref.k_bound), order=m,
                total_entries=total_entries, cost=dict(ref.cost))


def construct(h,gamma,n,initial,duration,tolerance,max_order=24,cap=512):
    if type(max_order) is not int or not 0<=max_order<=24: raise ValueError('order cap')
    duration,tolerance=rational(duration),rational(tolerance)
    if duration<=0 or tolerance<0: raise ValueError('duration/tolerance')
    initial=clean(initial,n,cap)
    ref=Reference(h,gamma,n,cap)
    layers=[ref.initial(initial)]
    probes=[]
    total_entries=len(layers[0])
    try:
        bounds=[rational(sum(map(abs,initial.values()),F(0))*(ref.k_bound*duration)**(k+1)/factorial(k+1)) for k in range(max_order+1)]
        if not any(b<=tolerance for b in bounds):raise BudgetExceeded('factorial preflight exceeds order budget')
        for k,bound in enumerate(bounds):
            probes.append(dict(order=k,bound=str(bound),entries=len(layers[-1])))
            if bound<=tolerance:
                return layers,dict(probes=probes,work=dict(ref.cost),total_entries=total_entries)
            if k<max_order:
                layer=ref.integrate(ref.apply_b(layers[-1]))
                total_entries+=len(layer)
                if total_entries>50000: raise BudgetExceeded('expansion total entry budget')
                layers.append(layer)
        raise BudgetExceeded('expansion order budget without certificate')
    except ValueError as exc:
        exc.construction_cost=dict(probes=probes,work=dict(ref.cost),total_entries=total_entries)
        raise


def evaluate_expansion(layers,n,duration,tolerance,cap=512):
    from .v10_exp_readout import evaluate_mode_map
    mapping={}
    for layer in layers:
        for (p,a,b,k),(cr,ci) in layer.items():mapping.setdefault(p,[]).append((a,b,k,cr,ci))
    mids,error,cost=evaluate_mode_map(mapping,duration,tolerance)
    pauli=eigen_to_pauli(mids,n,cap)
    # True Pm(T) is Hermitian (established by check_expansion). Symmetrizing
    # numerical A cannot increase operator-norm error to Hermitian Pm(T).
    out=clean({p:c[0] for p,c in pauli.items() if c[0]},n,cap)
    cost['output_pauli_terms']=len(out)
    cost['basis_conversion_paths']=sum(2**sum(c in 'PM' for c in p) for p in mids)
    return out,error,cost
