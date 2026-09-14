"""Inspect bounded cubic CAR motifs against an exact restricted dual."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import canonical, mono, product, encode
from experiments.marginal_hunt_car import adj
from research.certificate_scaling.commutator_dual_witness import moment_decode, gram, evaluate, check


def cubic_words(support, charge=-1):
    if charge == -3:
        return [tuple((0,i) for i in c) for c in combinations(support,3)]
    if charge != -1:
        raise ValueError('Negative charge -1 or -3 required')
    return [((0,i),) for i in support]+[
        ((1,i),(0,j),(0,k)) for i in support for j,k in combinations(support,2)]


def scan(witness_path,out):
    started=time.monotonic()
    raw=witness_path.read_bytes();witness=json.loads(raw)
    accepted=check(witness)
    checked=time.monotonic()
    m=witness['modes'];y=moment_decode(witness['moments'],m)
    import numpy as np
    families=[]
    for width in (3,4):
        for charge in (-1,-3):
            stats=[];pair_products=0
            for support in combinations(range(m),width):
                words=cubic_words(support,charge)
                polys=[mono(w) for w in words]
                g=gram(polys,y);pair_products+=len(words)*(len(words)+1)
                ev,u=np.linalg.eigh(np.array(g,dtype=float))
                v=[F(int(round(float(x)*10**6)),10**6) for x in u[:,0]]
                exact=sum(v[i]*g[i][j]*v[j] for i in range(len(words)) for j in range(len(words)))
                if exact < 0:
                    p={w:c for w,c in zip(words,v) if c}
                    square=product(canonical(adj(p)),p)
                    if evaluate(square,y)!=exact:raise AssertionError('Rayleigh/CAR mismatch')
                    stats.append({'support':support,'numeric_minimum':float(ev[0]),
                                  'exact_negative_square':str(exact),'operator':encode(p),
                                  'nonzeros':len(p),'minimum_diagonal':str(min(g[i][i] for i in range(len(words))))})
            stats.sort(key=lambda x:x['numeric_minimum'])
            families.append({'width':width,'charge':charge,'clusters_scanned':len(list(combinations(range(m),width))),
                             'word_pair_products':pair_products,'negative_clusters':len(stats),'best':stats[:8]})
    result={'witness_sha256':hashlib.sha256(raw).hexdigest(),'dual_replay':accepted,
            'families':families,'dual_replay_seconds':checked-started,
            'scan_seconds':time.monotonic()-checked,'wall_seconds':time.monotonic()-started,
            'scope':'Exact negative directions exclude this dual from a larger cone; they do not by themselves improve a primal lower bound.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in row.items() if k!='best'}|{'best_value':row['best'][0]['numeric_minimum'] if row['best'] else None} for row in families]),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--witness',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();scan(a.witness,a.out)
