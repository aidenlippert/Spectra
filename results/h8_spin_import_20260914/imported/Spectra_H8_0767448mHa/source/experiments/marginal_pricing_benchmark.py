"""Compare complete and stabilizer-split moment pricing on fixed duals."""
import json
from pathlib import Path
import time

import numpy as np

from experiments.marginal_energy_selector import prepare, price_candidates, restricted_problem


def run(modes=10):
    started=time.monotonic();data=prepare(modes,pricing='stabilizer')
    build=time.monotonic()-started;full=dict(data,pricing='full')
    dual=restricted_problem(data,[],8)()['dual']
    rng=np.random.default_rng(692)
    duals=[dual]+[rng.normal(size=len(dual)) for _ in range(19)]
    eigenvalue_error=0.;column_error=0.
    for y in duals:
        a,na=price_candidates(full,y,1000);b,nb=price_candidates(data,y,1000)
        if na!=nb or {p['candidate'] for p in a}!={p['candidate'] for p in b}:
            raise AssertionError('Pricing dictionaries disagree')
        for p in b:
            i=p['candidate'];v=p['vector']
            eigenvalue_error=max(eigenvalue_error,abs(p['violation']-next(q['violation'] for q in a if q['candidate']==i)))
            column_error=max(column_error,float(np.max(abs(data['candidate_maps'][i]@np.outer(v,v).ravel()-p['column']))))
    if eigenvalue_error>1e-9 or column_error>1e-9:raise AssertionError('Pricing equivalence failed')
    timings={'full':[],'stabilizer':[]}
    for k in range(5):
        for label in (['full','stabilizer'] if k%2==0 else ['stabilizer','full']):
            started=time.perf_counter()
            for y in duals:price_candidates(full if label=='full' else data,y,8)
            timings[label].append((time.perf_counter()-started)/len(duals))
    receipt={'modes':modes,'build_seconds':build,'full_candidate_blocks':len(data['candidates']),
             'split_candidate_blocks':len(data['pricing_blocks']),
             'largest_full_matrix':max(len(b['words']) for b in data['candidates']),
             'largest_split_matrix':max(b['transform'].shape[1] for b in data['pricing_blocks']),
             'maximum_minimum_eigenvalue_difference':eigenvalue_error,'maximum_lifted_column_difference':column_error,
             'pricing_seconds_per_call':timings,
             'scope':'Alternating batches on 20 fixed dual vectors; includes contraction, diagonalization, and column construction. Full maps and permutation enumeration still used during preparation.'}
    out=Path(__file__).resolve().parents[1]/'results/marginal_reduced_pricing';out.mkdir(exist_ok=True)
    (out/f'm{modes}_benchmark.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--modes',type=int,default=10)
    run(parser.parse_args().modes)
