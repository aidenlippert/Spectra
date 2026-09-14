"""Independent exact 70x70 fixed-N oracle for the exported rational H4 CAR fixture."""
from fractions import Fraction as F
from itertools import combinations
import json, math
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[2].parents[0]
FIX = ROOT / 'results/certificate_scaling/active_space_ladder/h4/fixture.json'
LEGACY = ROOT / 'results/marginal_molecule/h4_rectangle_sto3g.json'

def matrix(data, particles=4):
    m = data['modes']; states=[sum(1<<i for i in c) for c in combinations(range(m),particles)]
    pos={s:i for i,s in enumerate(states)}; n=len(states); A=[[F(0) for _ in range(n)] for _ in range(n)]
    for item in data['hamiltonian']:
        word=[tuple(x) for x in item['word']]; c=F(item['coefficient']);
        for j,s0 in enumerate(states):
            s=s0; sign=1
            for creation,i in reversed(word):
                sign *= -1 if ((s & ((1<<i)-1)).bit_count() & 1) else 1
                if creation:
                    if s>>i & 1: break
                    s |= 1<<i
                else:
                    if not (s>>i & 1): break
                    s &= ~(1<<i)
            else: A[pos[s]][j] += c*sign
    assert all(A[i][j] == A[j][i] for i in range(n) for j in range(n))
    return A, states

def ldl_psd(A):
    """Exact fraction-free-enough LDL test; all pivots must be nonnegative."""
    n=len(A); L=[[F(0)]*n for _ in range(n)]; D=[F(0)]*n
    for i in range(n):
        for j in range(i):
            numer=A[i][j]-sum(L[i][k]*D[k]*L[j][k] for k in range(j))
            if D[j] == 0:
                if numer != 0: return False, F(-1), j
                L[i][j]=F(0)
            else: L[i][j]=numer/D[j]
        D[i]=A[i][i]-sum(L[i][k]*L[i][k]*D[k] for k in range(i))
        if D[i] < 0: return False, D[i], i
        L[i][i]=F(1)
    return True, min(D,default=F(0)), None

def rayleigh(A, v):
    den=sum(x*x for x in v); return sum(v[i]*A[i][j]*v[j] for i in range(len(v)) for j in range(len(v)))/den

def main():
    raw=FIX.read_bytes(); source_sha=hashlib.sha256(raw).hexdigest()
    assert source_sha == '8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120'
    data=json.loads(raw); assert data['modes']==8 and data['particles']==4
    A,states=matrix(data); n=len(A); assert n == 70
    import numpy as np
    vals,vecs=np.linalg.eigh(np.array([[float(x) for x in row] for row in A])); e=float(vals[0])
    # Search downward rational shifts until exact PSD succeeds.
    scale=10**8; lower=F(math.floor(e*scale)-1,scale)
    while True:
        B=[[A[i][j]-(lower if i==j else 0) for j in range(n)] for i in range(n)]
        ok,pivot,idx=ldl_psd(B)
        if ok: break
        lower-=F(1,scale)
    v=[F(int(round(float(x)*10**9)),10**9) for x in vecs[:,0]]
    upper=rayleigh(A,v)
    assert lower <= upper
    # Regression: zero pivot with coupling is indefinite and must fail.
    assert not ldl_psd([[F(0),F(1)],[F(1),F(0)]])[0]
    # Negative controls: a corrupted matrix and a same-size wrong fixture must fail.
    bad=[row[:] for row in A]; bad[0][0] -= F(1,100)
    assert not ldl_psd([[bad[i][j]-(lower if i==j else 0) for j in range(n)] for i in range(n)])[0]
    wrong=json.loads(LEGACY.read_text()); _, wrong_states=matrix(wrong)
    assert len(wrong_states) == n and wrong['modes']==data['modes'] and wrong['particles']==data['particles']
    assert hashlib.sha256(LEGACY.read_bytes()).hexdigest() != source_sha
    out={'dimension':n,'sector_particles':4,'source_sha256':source_sha,'legacy_stress_sha256':hashlib.sha256(LEGACY.read_bytes()).hexdigest(),'numerical_ground':e,'lower':str(lower),'upper':str(upper),'width':str(upper-lower),'lower_float':float(lower),'upper_float':float(upper),'width_float':float(upper-lower),'min_ldl_pivot':str(pivot),'work':'O(C(8,4)^3) exact rational LDL; baseline oracle, not scalable discovery','negative_controls':'corrupted diagonal rejected; same-dimension legacy fixture rejected by content hash'}
    out_path=ROOT/'results/all_angles_20260913/frontier/h4_exact_oracle.json'; out_path.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:out[k] for k in ('dimension','lower_float','upper_float','width_float','negative_controls')},indent=2)); return out

if __name__=='__main__': main()
