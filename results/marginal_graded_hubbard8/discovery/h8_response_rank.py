"""Exact rank obstruction for the H8 valence response budget.

The calculation is deliberately independent of response construction: it builds
the 70 (4 up, 4 down) determinants, forms the exact leakage Gram matrix, and
checks principal minors of the scalar Schur matrix.
"""
from fractions import Fraction as F
from itertools import combinations
import argparse, json, hashlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_schur_transfer import ldl_pivots

ROOT = Path(__file__).resolve().parents[1]

def valence_states():
    out=[]
    for up in combinations(range(8),4):
        # One electron per site: the down set is the complement of up.
        out.append(sum(1 << (2*i) for i in up) |
                   sum(1 << (2*i+1) for i in range(8) if i not in up))
    return out

def gram_leakage(oracle, states):
    retained=set(states); cols=[]
    for s in states:
        cols.append({t:v for t,v in oracle.action(s).items() if t not in retained})
    g=[]
    for a in cols:
        g.append([sum(x*b.get(t,F(0)) for t,x in a.items()) for b in cols])
    return g

def find_principal(s, want=33):
    import numpy as np
    a=np.asarray(s,dtype=float)
    vals,vecs=np.linalg.eigh(a)
    neg=int(np.sum(vals < -1e-8))
    # A deterministic greedy coordinate search, followed by exact verification.
    order=list(np.argsort(np.diag(a))[::-1])
    for seed in (order, list(np.argsort(vals))):
        chosen=[]
        for i in seed:
            trial=chosen+[int(i)]
            if np.linalg.eigvalsh(a[np.ix_(trial,trial)])[ -1] < -1e-8:
                chosen=trial
            if len(chosen)>=want: return chosen,neg,vals
    return chosen,neg,vals

def main(out, taus):
    cert=json.loads((ROOT/'hamiltonian.json').read_text())
    cert['modes']=16; cert['particles']=8
    oracle=DeterminantOracle(cert); p=valence_states(); g=gram_leakage(oracle,p)
    gamma=F(-413,100); result={'schema':'h8-response-rank-v1','gamma':str(gamma),
      'retained_states':p,'retained_dimension':len(p),'gram':[[str(x) for x in r] for r in g],
      'scope':'Exact CAR action and rational leakage Gram for fixed-up4/down4 valence determinants.'}
    result['hamiltonian_canonical_sha256']=hashlib.sha256(json.dumps(cert,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    rows=[]
    for tau_text in taus:
        tau=F(tau_text); den=gamma-tau
        s=[[-tau*(i==j)-g[i][j]/den for j in range(len(p))] for i in range(len(p))]
        inds,neg,vals=find_principal(s)
        exact=ldl_pivots([[-s[i][j] for j in inds] for i in inds]) if len(inds)>=33 else None
        row={'tau':str(tau),'negative_inertia':neg,'principal_indices':inds,
             'principal_dimension':len(inds),'exact_ldl_positive':exact is not None,
             'ldl_pivots':[str(x) for x in exact] if exact else []}
        rows.append(row)
    result['tests']=rows
    result['note']='A negative 33-dimensional restriction of S means every PSD response correction certifying this tau needs rank at least 33; this is a response-budget obstruction, not a physical gap claim.'
    out=Path(out); out.mkdir(parents=True,exist_ok=True); (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'summary.json').write_text(json.dumps({'tests':rows,'action_cache':len(oracle.cache),'referenced_determinants':oracle.referenced_state_count()},indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--out',default=str(ROOT/'response_rank')); ap.add_argument('--tau',action='append',default=['-4.2','-4.25','-4.3']); a=ap.parse_args(); main(a.out,a.tau)
