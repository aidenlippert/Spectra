"""Bounded exact recognizer for matching-dressed quadratic Hamiltonians."""
from collections import deque
from fractions import Fraction

from experiments.marginal_hunt_car import add, mul, mono, scale
from experiments.marginal_symbolic import canonical, hermitian


def _prod(*ps):
    out = {(): Fraction(1)}
    for p in ps:
        out = mul(out, p)
    return out


def _num(i):
    return mono(((1, i), (0, i)))


def dress_polynomial(poly, matching, m):
    """Conjugate a normal-ordered CAR polynomial by a bounded matching."""
    if type(m) is not int or m < 0: raise ValueError("invalid mode count")
    if not isinstance(poly, dict): raise ValueError("invalid polynomial")
    for word, coeff in poly.items():
        if type(coeff) not in (int, Fraction) or type(word) is not tuple or any(
            type(t) is not tuple or len(t)!=2 or type(t[0]) is not int or t[0] not in (0,1) or
            type(t[1]) is not int or not 0<=t[1]<m for t in word):
            raise ValueError("invalid exact polynomial")
    mate = {}
    for pair in matching:
        if type(pair) is not tuple or len(pair) != 2: raise ValueError("invalid matching")
        i,j=pair
        if type(i) is not int or type(j) is not int or not (0<=i<m and 0<=j<m) or i==j or i in mate or j in mate:
            raise ValueError("invalid matching")
        mate[i]=j; mate[j]=i
    controls = {i: mate[i] for i in mate}
    out = {}
    for word, coeff in poly.items():
        factors = []
        for creation, i in word:
            if not 0 <= i < m:
                raise ValueError("mode outside range")
            p = {(): Fraction(1)}
            if i in controls: p = add(p, scale(_num(controls[i]), -2))
            factors.append(_prod(p, mono(((1, i),))) if creation else
                          _prod(mono(((0, i),)), p))
        out = add(out, scale(_prod(*factors), coeff))
    return out


def _canonical(poly):
    return {w: c for w, c in poly.items() if c}


def recognize(h, m):
    if type(m) is not int or m < 4: raise ValueError("m must be at least 4")
    if not isinstance(h, dict): return {"accepted": False, "reason": "invalid polynomial"}
    for w,c in h.items():
        if type(c) is bool or not isinstance(c,(int,Fraction)): return {"accepted":False,"reason":"non-rational coefficient"}
        if type(w) is not tuple or len(w) not in (0,2,4) or any(type(x) is not tuple or len(x)!=2 or type(x[0]) is not int or x[0] not in (0,1) or type(x[1]) is not int or not 0<=x[1]<m for x in w):
            return {"accepted": False, "reason": "invalid polynomial word"}
    if any(len(w) and sum(x[0] for x in w)*2 != len(w) for w in h): return {"accepted":False,"reason":"non-number-conserving"}
    if any(canonical({w: Fraction(1)}) != {w: Fraction(1)} for w in h):
        return {"accepted":False,"reason":"noncanonical input word"}
    h = {w:Fraction(c) for w,c in h.items() if c}
    if not hermitian(h): return {"accepted":False,"reason":"non-Hermitian"}
    q, quart = {}, {}
    for w, c in h.items():
        c = Fraction(c)
        if len(w) in (0, 2): q[w] = q.get(w, 0) + c
        else: quart[w] = quart.get(w, 0) + c
    for w in q:
        if not w: continue
        if w[0][0] != 1 or w[1][0] != 0:
            return {"accepted": False, "reason": "unsupported quadratic term"}
    if not quart:
        return {"accepted": True, "matching": (), "recovered_quadratic": q,
                "graph_work": 0, "bfs_work": 0, "exact_regeneration_match": True}
    edges = {(w[0][1], w[1][1]) for w,c in q.items() if c and len(w)==2 and w[0][1]!=w[1][1]}
    adj = {i:set() for i in range(m)}
    for i,j in edges: adj[i].add(j); adj[j].add(i)
    r = {}
    for w,c in quart.items():
        cs = [i for a,i in w if a]; an = [i for a,i in w if not a]
        common = set(cs)&set(an)
        if len(cs)!=2 or len(an)!=2 or len(common)!=1: return {"accepted":False,"reason":"unsupported quartic term"}
        k=next(iter(common)); i=next(x for x in cs if x!=k); j=next(x for x in an if x!=k)
        if i==j or not q.get(((1,i),(0,j)), 0): return {"accepted":False,"reason":"quartic without hopping"}
        key=(i,j)
        if key in r: return {"accepted":False,"reason":"multiple controls"}
        # The unique control coefficient is -2*h_ij times its CAR sign.
        expected = scale(mul(mono(((1,i),(0,j))), _num(k)), -2*q[((1,i),(0,j))])
        if expected != {w:c}: return {"accepted":False,"reason":"invalid controlled hopping coefficient"}
        r[key]=k
    # Infer one XOR bit per ordered hop, then solve column offsets by BFS.
    graph_work=0; bfs_work=0
    base = {}
    for k in range(m):
        seen={}; root=next((i for i in range(m) if i!=k), None); seen[root]=0; dq=deque([root])
        while dq:
            i=dq.popleft(); bfs_work+=1
            for j in adj[i]:
                if j==k: continue
                graph_work+=1; val=seen[i]^(int(r.get((i,j))==k))
                if j in seen and seen[j]!=val: return {"accepted":False,"reason":"inconsistent BFS"}
                if j not in seen: seen[j]=val; dq.append(j)
        if len(seen)<m-1: return {"accepted":False,"reason":"vertex-deletion-disconnected"}
        for i,v in seen.items(): base[i,k]=v
    # r_ik xor r_ki is the symmetric constraint; solve offsets.
    g={0:0}
    for k in range(1,m): g[k]=base[0,k]^base[k,0]
    if any((g[i]^g[k]) != (base[i,k]^base[k,i]) for i in range(m) for k in range(i+1,m)):
        return {"accepted":False,"reason":"inconsistent symmetric offsets"}
    candidates=[]
    for flip in (0,1):
        S=tuple(tuple(sorted(j for j in range(m) if j!=i and base.get((i,j),0)^g.get(j,0)^flip)) for i in range(m))
        if all((j in S[i]) == (i in S[j]) for i in range(m) for j in range(m)) and all(len(x)<=1 for x in S): candidates.append(S)
    for S in candidates:
        matching=tuple((i,S[i][0]) for i in range(m) if S[i] and i<S[i][0])
        if len({x for p in matching for x in p}) != 2*len(matching): continue
        pred=dress_polynomial(q, matching, m)
        if _canonical(pred)==_canonical(h):
            return {"accepted":True,"matching":matching,"recovered_quadratic":q,"graph_work":graph_work,"bfs_work":bfs_work,"exact_regeneration_match":True}
    return {"accepted":False,"reason":"no matching candidate"}
