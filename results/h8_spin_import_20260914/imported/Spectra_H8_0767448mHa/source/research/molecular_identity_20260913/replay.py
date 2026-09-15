"""Standard-library-only replay of exact local obstruction directions."""
from fractions import Fraction as F
from itertools import product as bit_patterns
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import canonical,decode,mono,add,scale,product
from experiments.marginal_hunt_car import adj
from research.certificate_scaling.commutator_dual_witness import check,moment_decode,evaluate


def occupation_projector(support,bits):
    if len(support)!=len(bits) or len(set(support))!=len(support) or any(type(i) is not int or i<0 for i in support) or any(type(b) is not int or b not in (0,1) for b in bits):
        raise ValueError('Distinct orbital labels and binary occupancies required')
    operator=canonical(mono(tuple((1-b,i) for i,b in zip(support,bits))))
    projector=mono(())
    for i,b in zip(support,bits):
        ni=mono(((1,i),(0,i)))
        projector=product(projector,ni if b else add(mono(()),scale(ni,-1)))
    if product(canonical(adj(operator)),operator)!=projector:
        raise AssertionError('Occupation square identity failed')
    return operator,projector


def replay(scan,witness_raw):
    started=time.monotonic()
    if hashlib.sha256(witness_raw).hexdigest()!=scan['witness_sha256']:
        raise ValueError('Source witness hash mismatch')
    witness=json.loads(witness_raw);dual=check(witness)
    m=witness['modes'];y=moment_decode(witness['moments'],m);count=0;minimum=F(0)
    for family in scan['families']:
        for example in family['best']:
            support=example['support']
            if len(support)!=family['width'] or len(set(support))!=len(support) or any(type(i) is not int or not 0<=i<m for i in support):
                raise ValueError('Invalid local support')
            p=decode(example['operator'],m,3)
            if not p or any(i not in support for w in p for _,i in w):raise ValueError('Operator leaves local support')
            if {sum(2*f-1 for f,_ in w) for w in p}!={family['charge']}:raise ValueError('Wrong operator charge')
            value=evaluate(product(canonical(adj(p)),p),y)
            if value!=F(example['exact_negative_square']) or value>=0:
                raise ValueError('Claimed exact negative direction failed')
            minimum=min(minimum,value);count+=1
    support=(2,4,5)
    identities=[occupation_projector(support,bits)[1] for bits in bit_patterns((0,1),repeat=3)]
    if add(*identities)!=mono(()):raise AssertionError('Occupation projectors do not sum to identity')
    _,p=occupation_projector(support,(1,1,1));value=evaluate(p,y)
    if value>=0:raise ValueError('Three-occupation obstruction is absent')
    return {'accepted_local_negative_directions':count,'minimum_exact_value':str(minimum),
            'negative_occupation_support':support,'negative_occupation_value':str(value),
            'eight_occupation_projector_identities_checked':True,
            'source_dual_objective':dual['dual_objective'],
            'source_witness_sha256':hashlib.sha256(witness_raw).hexdigest(),
            'many_body_space_enumerated':False,'replay_seconds':time.monotonic()-started,
            'scope':'Exact separation of this restricted dual by valid local CAR squares. No strict improvement of the optimized energy bound follows from separation alone.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scan',type=Path,required=True)
    p.add_argument('--witness',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=replay(json.loads(a.scan.read_text()),a.witness.read_bytes())
    a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
