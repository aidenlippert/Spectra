"""Bounded consistency diagnostic; discovery uses no higher moments."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse
import hashlib
import json
import time
import numpy as np

from experiments.marginal_symbolic import mono,canonical,add,scale
from research.joint_patterns_20260913.dual import integer_grams
from research.joint_patterns_20260913.core import anticommutator
from research.certificate_scaling.commutator_dual_witness import moment_decode,evaluate

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'results/response_consistency_20260913'


def run(out=OUT):
    out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();source=ROOT/'results/trace_pricing_20260913/full_dual/witness.json'
    witness=json.loads(source.read_text());y=moment_decode(witness['moments'],12)
    triples=list(combinations(range(12),3));polys=[canonical(mono(tuple((0,i) for i in t))) for t in triples]
    gram,den=integer_grams(polys,[y],'anticommutator');G=np.array([[float(F(x,den)) for x in row] for row in gram[0]])
    diagonal=min((G[i,i],t) for i,t in enumerate(triples));records=[];chosen=None;rounding=[]
    for spatial_count in (3,4,5,6):
        proposals=[]
        for subset in combinations(range(6),spatial_count):
            ids=[i for i,t in enumerate(triples) if all(j//2 in subset for j in t)]
            v,U=np.linalg.eigh(G[np.ix_(ids,ids)])
            proposals.append((float(v[0]),subset,ids,U[:,0]))
        best=min(proposals,key=lambda x:x[0]);value,subset,ids,v=best
        records.append({'spatial_orbitals':spatial_count,'subsets_checked':len(proposals),'minimum_eigenvalue':value,'best_support':subset})
        if value < -1e-9:
            for denominator in (10,100,1000,10**6):
                integers=np.rint(v*denominator).astype(int)
                C=add(*(scale(polys[i],F(int(z),denominator)) for i,z in zip(ids,integers) if z))
                P=anticommutator(C,C)
                if any(len(w)>4 or sum(2*c-1 for c,_ in w) for w in P):raise AssertionError('Outside the specified moment domain')
                moment=evaluate(P,y)
                rounding.append({'denominator':denominator,'exact_moment':str(moment),'accepted':moment<0,
                    'triples':[list(triples[i]) for i,z in zip(ids,integers) if z],
                    'coefficients':[int(z) for z in integers if z]})
                if moment<0:
                    chosen={'kind':'three_removal_separator_v1','parent_witness_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                        'modes':12,'support_spatial':list(subset),'triples':[list(triples[i]) for i,z in zip(ids,integers) if z],
                        'coefficients':[int(z) for z in integers if z],'denominator':denominator,
                        'exact_moment':str(moment),'moment_float':float(moment),'quartic_terms':len(P),
                        'maximum_degree':max(map(len,P)), 'scope':'Positive anticommutator with exactly cancelled sixth degree; only the existing degree-four functional was evaluated.'}
                    break
        if chosen:break
    result={'diagonal_minimum':{'value':diagonal[0],'triple':diagonal[1]},'search':records,'separator':chosen,'rounding_attempts':rounding,
        'wall_seconds':time.monotonic()-start,'many_body_states_enumerated':0}
    (out/'t1_probe.json').write_text(json.dumps(result,indent=2)+'\n')
    if chosen:(out/'separator.json').write_text(json.dumps(chosen,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=OUT)
    run(parser.parse_args().out)
