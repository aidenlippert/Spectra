"""Exact verifier for rational SOS orbit certificates."""
from fractions import Fraction as F
from experiments.marginal_symbolic import decode, product, add, scale, mono, adj, hermitian, number_shift

def permuted_canonical_word(word, mapping):
    """Signed image of a canonical normal-ordered monomial."""
    left=[mapping[i] for c,i in word if c]
    right=[mapping[i] for c,i in word if not c]
    parity=sum(a>b for seq in (left,right) for j,a in enumerate(seq) for b in seq[j+1:])
    image=tuple((1,i) for i in sorted(left))+tuple((0,i) for i in sorted(right))
    return image,(-1 if parity%2 else 1)


def average_canonical_polynomial(poly, permutations):
    # Count integer signed images before doing large rational arithmetic.
    result={}
    for word,coefficient in poly.items():
        counts={}
        for mapping in permutations:
            image,sign=permuted_canonical_word(word,mapping)
            counts[image]=counts.get(image,0)+sign
        for image,count in counts.items():
            if count:result[image]=result.get(image,0)+coefficient*F(count,len(permutations))
    return {w:c for w,c in result.items() if c}

def verify(c):
    m,n=c.get('modes'),c.get('particles'); d=c.get('operator_degree',3)
    if type(m) is not int or type(n) is not int or not 1<=m or not 0<=n<=m: raise ValueError('Invalid sector')
    if type(d) is not int or d not in (3,4): raise ValueError('Invalid degree')
    h=decode(c['hamiltonian'],m,4); x=decode(c['number_multiplier'],m,2*d-2)
    if not hermitian(h) or not hermitian(x): raise ValueError('H and X must be Hermitian')
    if any(sum(2*cr-1 for cr,_ in w) for w in h) or any(sum(2*cr-1 for cr,_ in w) for w in x): raise ValueError('Nonconserving polynomial')
    if type(c.get('orbit_squares')) is not list or type(c.get('b')) is not str: raise ValueError('Invalid certificate fields')
    perms=c.get('permutations')
    if not isinstance(perms,list) or not perms: raise ValueError('Missing permutations')
    P=[]
    for q in perms:
        if type(q) is not list or len(q)!=m or any(type(i) is not int for i in q) or sorted(q)!=list(range(m)): raise ValueError('Invalid permutation')
        P.append(q)
    def square_sum(items):
        if type(items) is not list:raise ValueError('Invalid square list')
        result={}
        for item in items:
            if not isinstance(item,dict): raise ValueError('Invalid square')
            p=decode(item.get('polynomial',[]),m,d)
            if not p: raise ValueError('Empty polynomial')
            charges={sum(2*cr-1 for cr,_ in w) for w in p}
            if len(charges)!=1: raise ValueError('Mixed-charge polynomial')
            if type(item.get('weight')) is not str: raise ValueError('Weight must be rational string')
            weight=F(item['weight'])
            if weight<0: raise ValueError('Negative weight')
            result=add(result,scale(product(adj(p),p),weight))
        return result
    # Positive orbit averages remain valid when H has less symmetry. Local
    # additions must stay outside that average to represent asymmetric terms.
    sos=add(average_canonical_polynomial(square_sum(c['orbit_squares']),P),
            square_sum(c.get('direct_squares',[])))
    b=F(c['b']); residual=add(h,mono((),-b),scale(sos,-1),scale(product(number_shift(m,n),x),-1))
    if not hermitian(residual): raise ValueError('Non-Hermitian residual')
    eta=sum(abs(v) for v in residual.values())
    return {'lower':str(b-eta),'lower_float':float(b-eta),'residual_l1':str(eta),'residual_l1_float':float(eta),
            'residual_constant':str(residual.get((),F(0))),
            'residual_nonconstant_l1':str(eta-abs(residual.get((),F(0)))),
            'orbit_squares':len(c['orbit_squares']),'direct_squares':len(c.get('direct_squares',[])),
            'permutations':len(P),'residual_terms':len(residual)}


def normalize_scalar_residual(certificate):
    """Absorb the residual's scalar into b and recheck the new certificate."""
    before=verify(certificate)
    normalized=dict(certificate)
    normalized['b']=str(F(certificate['b'])+F(before['residual_constant']))
    return normalized,verify(normalized)

if __name__=='__main__':
    import json,sys
    print(json.dumps(verify(json.load(open(sys.argv[1]))),indent=2))
