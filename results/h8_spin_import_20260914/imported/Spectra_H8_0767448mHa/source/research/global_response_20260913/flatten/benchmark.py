"""Bounded execution comparison on the largest charged H6 conserved block."""
import json, time
from pathlib import Path
import numpy as np
from research.compact_response_20260913 import program, closure
from .flatten import nested_action, flattened_action, schedule_counts

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'results/global_response_20260913/flatten'

def main():
    data,_,_=program.load_case('h6'); blocks,den,prep=closure.blocks(data)
    block=max(blocks,key=lambda x:len(x['H'])); n=len(block['H']); H=np.asarray(block['H'],dtype=float)/den
    first=json.loads((program.OUT/'h6_program.json').read_text())
    from fractions import Fraction
    b=float(Fraction(first['target_Ha']))
    H-=b*np.eye(n)
    # Physical occupation partitions from the actual charged H6 states.
    states=block['states']; m=data['modes']
    v=np.zeros(n); v[states.index(63) if 63 in states else 0]=1.0
    q1mask=[((s>>(m-2))&3)==3 for s in states]
    r2mask=[(not q1mask[i]) and (((s>>(m-4))&3)==3) for i,s in enumerate(states)]
    p2mask=[(not q1mask[i]) and (not r2mask[i]) for i in range(n)]
    q1=sum(q1mask); q2=sum(r2mask)
    def mat(x): return H@np.asarray(x,dtype=float)
    def P1(x):
        y=np.asarray(x,dtype=float).copy(); y[np.asarray(q1mask)]=0.; return y.tolist()
    def Q1(x):
        y=np.zeros(n); y[np.asarray(q1mask)]=np.asarray(x)[np.asarray(q1mask)]; return y.tolist()
    def P2(x):
        y=np.zeros(n); y[np.asarray(p2mask)]=np.asarray(x)[np.asarray(p2mask)]; return y.tolist()
    def Q2(x):
        y=np.zeros(n); y[np.asarray(r2mask)]=np.asarray(x)[np.asarray(r2mask)]; return y.tolist()
    r1=first
    joint=json.loads((ROOT/'results/composable_response_20260913/h6_joint_response.json').read_text())
    r2=joint['second_response']
    from fractions import Fraction
    eta=float(Fraction(program.check(data, json.loads((program.ROOT/'results/molecular_collective_20260913/campaign/h6/rank_10/tail.json').read_text()), r1)['exact_program_residual_penalty_Ha']))
    from research.composable_response_20260913 import recursive
    def K1(x): return recursive.retained_action(mat,P1,Q1,x,r1)
    eta2=float(Fraction(r2['residual_penalty_Ha']))
    def D2(x):
        y=np.asarray(K1(x),dtype=float); y-=eta*np.asarray(P1(x),dtype=float); return y.tolist()
    t=time.perf_counter(); ref=recursive.retained_action(D2,P2,Q2,v.tolist(),r2); ref=np.asarray(ref)-eta2*np.asarray(P2(v.tolist())); nested_s=time.perf_counter()-t
    from .flatten import MemoAction
    memo=MemoAction(mat)
    def K1f(x): return recursive.retained_action(memo,P1,Q1,x,r1)
    def D2f(x):
        y=np.asarray(K1f(x),dtype=float); y-=eta*np.asarray(P1(x),dtype=float); return y.tolist()
    t=time.perf_counter(); flat=recursive.retained_action(D2f,P2,Q2,v.tolist(),r2); flat=np.asarray(flat)-eta2*np.asarray(P2(v.tolist())); flat_s=time.perf_counter()-t
    stats={'unique_H_rhs':memo.calls,'H_rhs':memo.h_rhs,'cache_entries':len(memo.cache)}
    err=float(np.max(np.abs(np.asarray(ref)-np.asarray(flat))))
    result={'fixture':'molecular_collective_20260913/campaign/h6','block_key':block['key'],
      'block_dimension':n,'charged_block_preparation':prep,'q1_dimension':q1,'q2_dimension':q2,
      'nested_seconds':nested_s,'flattened_seconds':flat_s,'speedup_nested_over_flattened':nested_s/flat_s,
      'max_abs_difference':err,'formal_counts':schedule_counts(r1['order'],r2['order']),
      'flatten_counts':stats,'storage_vectors_cached':stats['cache_entries'],
      'precision':'float64 diagnostic only','many_body_states_constructed':len(block['states']),
      'eta1_shift_Ha':eta,'first_order':r1['order'],'second_order':r2['order'],
      'scope':'The block is an enumerated H6 diagnostic; this benchmark makes no enumeration-free claim.'}
    OUT.mkdir(parents=True,exist_ok=True); (OUT/'benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
