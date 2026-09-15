"""Hamiltonian-only real orbital commutant discovery and exact transport proof.

Numerical rank/eigenvectors are proposals. Accepted certificates replay an
exact rational orthogonal map, transformed Hamiltonian, and inner CAR proof.
This recognizes hidden members of a known solvable family, not general chemistry.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json, sys, time
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import canonical, mono, add, scale, encode, decode, hermitian
from research.certificate_scaling.fermionic_ratio_chain import hamiltonian, compile_h, structured_replay


def generators(m):
    return [(i,j) for i in range(m) for j in range(i,m)]


def commutator_column(v, i, j, antisymmetric=False):
    """[Q(Eii),V] or [Q(Eij+Eji),V], by the CAR derivation rule."""
    out = {}
    for w, coefficient in v.items():
        for k, (cr, p) in enumerate(w):
            if p != i and p != j:
                continue
            q = j if p == i else i
            image = w[:k] + ((cr,q),) + w[k+1:]
            for z, sign in canonical(mono(image)).items():
                weight = (1 if p==j else -1) if antisymmetric else (1 if cr else -1)
                out[z] = out.get(z,F(0)) + coefficient * sign * weight
    return {w:c for w,c in out.items() if c}


def commutator_map(h, m, antisymmetric=False):
    v = {w:c for w,c in h.items() if len(w)==4}
    basis=list(combinations(range(m),2)) if antisymmetric else generators(m)
    cols = [commutator_column(v,i,j,antisymmetric) for i,j in basis]
    words = sorted(set().union(*(set(c) for c in cols)))
    return cols, words


def modular_rank(rows, prime=1000003):
    """Rank over a finite field is an exact lower bound on rational rank."""
    a = [[(x.numerator % prime)*pow(x.denominator % prime,-1,prime)%prime for x in row] for row in rows]
    rank = 0
    for col in range(len(a[0]) if a else 0):
        pivot = next((j for j in range(rank,len(a)) if a[j][col]), None)
        if pivot is None: continue
        a[rank],a[pivot] = a[pivot],a[rank]
        inv = pow(a[rank][col],-1,prime)
        a[rank] = [v*inv%prime for v in a[rank]]
        for j in range(rank+1,len(a)):
            fac = a[j][col]
            if fac: a[j] = [(x-fac*y)%prime for x,y in zip(a[j],a[rank])]
        rank += 1
        if rank == len(a): break
    return rank


def detect(h, m, tolerance=1e-9):
    import numpy as np
    from scipy.linalg import qr
    start = time.monotonic()
    cols,words = commutator_map(h,m)
    index = {w:i for i,w in enumerate(words)}
    c = np.zeros((len(words),len(cols)))
    for j,col in enumerate(cols):
        for w,value in col.items(): c[index[w],j] = float(value)
    # Economy SVD needs at least as many rows as columns to include all null vectors.
    if len(words)<len(cols): c = np.vstack((c,np.zeros((len(cols)-len(words),len(cols)))))
    _,s,vh = np.linalg.svd(c,full_matrices=False)
    threshold = tolerance*max(float(s[0]) if len(s) else 0,1.)
    rank = int(np.count_nonzero(s>threshold))
    null = vh[rank:]
    # Select a small independent row minor numerically, then prove its rank mod p.
    _,_,pivots = qr(c.T,mode='economic',pivoting=True)
    selected = [words[int(i)] for i in pivots[:rank] if i<len(words)]
    exact_rank = modular_rank([[col.get(w,F(0)) for col in cols] for w in selected])
    mats = []
    for vector in null:
        a = np.zeros((m,m))
        for value,(i,j) in zip(vector,generators(m)): a[i,j]=a[j,i]=value
        mats.append(a)
    comm = max((float(np.linalg.norm(a@b-b@a)) for i,a in enumerate(mats) for b in mats[i+1:]),default=0.)
    rec = {'modes':m,'unknowns':len(cols),'equation_rows':len(words),'map_nonzeros':sum(map(len,cols)),
           'numeric_nullity':len(null),'exact_rank_lower_bound':exact_rank,
           'exact_nullity_upper_bound':len(cols)-exact_rank,'rank_minor_prime':1000003,
           'rank_minor_words':[[list(letter) for letter in w] for w in selected],
           'singular_values':s.tolist(),'nullspace_commutator_max_numeric':comm,
           'discovery_seconds':time.monotonic()-start,'many_body_states_enumerated':0,
           'rank_classification':'exact upper bound on nullity; numerical nullspace is not a proof'}
    return mats,rec


def replay_rank(h,rec):
    """Exact, stdlib replay of the reported rational rank lower bound."""
    m=rec['modes'];cols,_=commutator_map(h,m)
    words=[tuple(tuple(letter) for letter in w) for w in rec['rank_minor_words']]
    rank=modular_rank([[col.get(w,F(0)) for col in cols] for w in words],rec['rank_minor_prime'])
    if rank!=rec['exact_rank_lower_bound']: raise ValueError('Rank witness mismatch')
    return {'exact_rank_lower_bound':rank,'exact_nullity_upper_bound':len(cols)-rank,
            'rules_out_real_density_basis':len(cols)-rank<m}


def transpose(a): return [list(row) for row in zip(*a)]


def matmul(a,b):
    out = [[F(0) for _ in range(len(b[0]))] for _ in a]
    for i,row in enumerate(a):
        for k,x in enumerate(row):
            if x:
                for j,y in enumerate(b[k]):
                    if y: out[i][j] += x*y
    return out


def orthogonal(u):
    m=len(u)
    return bool(m) and all(len(row)==m for row in u) and matmul(transpose(u),u)==[[F(i==j) for j in range(m)] for i in range(m)]


def rotate(h,u):
    """Substitute d_i = sum_p U[p,i] a_p, using exterior-square matrices."""
    m=len(u); pairs=list(combinations(range(m),2)); ix={p:i for i,p in enumerate(pairs)}
    one=[[F(0) for _ in range(m)] for _ in range(m)]
    two=[[F(0) for _ in pairs] for _ in pairs]
    for w,c in canonical(h).items():
        if not w: continue
        if len(w)==2 and [cr for cr,_ in w]==[1,0]: one[w[0][1]][w[1][1]]+=c
        elif len(w)==4 and [cr for cr,_ in w]==[1,1,0,0]:
            two[ix[(w[0][1],w[1][1])]][ix[(w[2][1],w[3][1])]]+=c
        else: raise ValueError('Only number-conserving one/two-body real Hamiltonians supported')
    exterior=[[u[p][i]*u[q][j]-u[p][j]*u[q][i] for i,j in pairs] for p,q in pairs]
    one=matmul(matmul(u,one),transpose(u))
    two=matmul(matmul(exterior,two),transpose(exterior))
    result={():h[()]} if h.get(()) else {}
    for i in range(m):
        for j in range(m):
            if one[i][j]: result[((1,i),(0,j))]=one[i][j]
    for a,(i,j) in enumerate(pairs):
        for b,(k,l) in enumerate(pairs):
            if two[a][b]: result[((1,i),(1,j),(0,k),(0,l))]=two[a][b]
    return result


def recognize_path(h,m):
    """Recover ordering and hopping gauge from a density-basis Hamiltonian."""
    graph=[set() for _ in range(m)]
    for w,c in h.items():
        if len(w)==4:
            if w[0][1]!=w[2][1] or w[1][1]!=w[3][1]: raise ValueError('Recovered interaction is not density-density')
            i,j=w[0][1],w[1][1];graph[i].add(j);graph[j].add(i)
    endpoints=[i for i in range(m) if len(graph[i])==1]
    if len(endpoints)!=2 or any(len(g) not in (1,2) for g in graph): raise ValueError('Interaction graph is not a path')
    order=[endpoints[0]]
    while len(order)<m:
        nxt=graph[order[-1]]-set(order)
        if len(nxt)!=1: raise ValueError('Disconnected interaction graph')
        order.append(nxt.pop())
    signs=[F(1)]
    for i,j in zip(order,order[1:]):
        hopping=h.get(((1,i),(0,j)),F(0))
        if abs(hopping)!=1: raise ValueError('Ratio family requires unit magnitude hopping')
        signs.append(-signs[-1]*hopping)
    p=[[F(0) for _ in range(m)] for _ in range(m)]
    for k,i in enumerate(order):p[i][k]=signs[k]
    recovered=rotate(h,transpose(p))
    return recovered,p


def discover(h,m,n,max_denominator=1000000):
    import numpy as np
    start=time.monotonic();h=canonical(h)
    if not hermitian(h): raise ValueError('Non-Hermitian input')
    mats,rec=detect(h,m)
    if len(mats)!=m or rec['nullspace_commutator_max_numeric']>1e-6:
        raise ValueError('No M-dimensional commuting real-symmetric kernel detected')
    rng=np.random.default_rng(7919)
    errors=[]
    for attempt in range(4):
        combination=sum((c*a for c,a in zip(rng.normal(size=m),mats)))
        eigenvalues,u=np.linalg.eigh(combination)
        rational=[[F(float(v)).limit_denominator(max_denominator) for v in row] for row in u]
        if not orthogonal(rational): errors.append('rational orthogonality failed');continue
        try:
            recovered,p=recognize_path(rotate(h,transpose(rational)),m)
            inner,_=compile_h(recovered,m,n)
        except ValueError as exc: errors.append(str(exc));continue
        rotation=matmul(rational,p)
        cert={'kind':'exact_orthogonal_ratio_transport_v1','modes':m,'particles':n,
              'hamiltonian':encode(h),'rotation':[[str(v) for v in row] for row in rotation],
              'inner_certificate':inner}
        receipt=replay(cert)
        rec.update(receipt);rec.update({'total_discovery_seconds':time.monotonic()-start,
                    'rational_reconstruction_attempts':attempt+1,'max_denominator':max_denominator,
                    'candidate_joint_eigengap':float(min(np.diff(eigenvalues))),
                    'source_rotation_used':False,'source_factors_used':False})
        return cert,rec
    raise ValueError('Basis proposal did not exactly recognize family: '+repr(errors))


def replay(cert):
    start=time.monotonic()
    if cert['kind']!='exact_orthogonal_ratio_transport_v1':raise ValueError('Unknown proof kind')
    m,n=cert['modes'],cert['particles'];inner=cert['inner_certificate']
    if (m,n)!=(inner['modes'],inner['particles']):raise ValueError('Sector mismatch')
    u=[[F(v) for v in row] for row in cert['rotation']]
    if len(u)!=m or not orthogonal(u):raise ValueError('Rotation is not exactly orthogonal')
    h=decode(cert['hamiltonian'],m,4)
    expected=rotate(decode(inner['hamiltonian'],m,4),u)
    if h!=expected:raise ValueError('Hamiltonian transport mismatch')
    rec=structured_replay(inner)
    rec.update({'exact_orthogonality':True,'exact_hamiltonian_transport':True,
                'transport_replay_seconds':time.monotonic()-start,
                'rotation_entries':m*m,'observed_hamiltonian_terms':len(h),
                'proof_scope':'rational real orthogonal transport of known ratio family'})
    return rec


def hidden_fixture(m):
    """Fixture construction only; the discovering function receives H, M, N."""
    def householder(v):
        norm=sum(x*x for x in v)
        return [[F(i==j)-F(2*v[i]*v[j],norm) for j in range(m)] for i in range(m)]
    u=matmul(householder([i%5+1 for i in range(m)]),householder([(3*i+2)%7+1 for i in range(m)]))
    assert orthogonal(u)
    return rotate(hamiltonian(m,F(4)),u)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',type=Path);p.add_argument('--fixture',type=Path);p.add_argument('--rank-replay',type=Path)
    p.add_argument('--sizes',type=int,nargs='+',default=[4,6,8]);p.add_argument('--out',type=Path,default=Path('results/certificate_scaling/hidden_density_basis'))
    a=p.parse_args()
    if a.verify:print(json.dumps(replay(json.loads(a.verify.read_text())),indent=2))
    elif a.fixture:
        fixture=json.loads(a.fixture.read_text());h=decode(fixture['hamiltonian'],fixture['modes'],4)
        if a.rank_replay:
            print(json.dumps(replay_rank(h,json.loads(a.rank_replay.read_text())),indent=2));sys.exit(0)
        _,rec=detect(h,fixture['modes']);rec['fixture']=str(a.fixture)
        rec['fixture_sha256']=hashlib.sha256(a.fixture.read_bytes()).hexdigest()
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
    else:
        a.out.mkdir(parents=True,exist_ok=True);rows=[]
        for m in a.sizes:
            start=time.monotonic();h=hidden_fixture(m);fixture_seconds=time.monotonic()-start
            cert,rec=discover(h,m,m//2);path=a.out/f'M{m}.json';path.write_text(json.dumps(cert,separators=(',',':'))+'\n')
            rec.update({'fixture_seconds':fixture_seconds,'certificate_bytes':path.stat().st_size,'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            rows.append(rec);(a.out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rec),flush=True)
