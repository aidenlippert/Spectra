"""Reconstruct every Chebyshev residual exactly using local integer contractions.

Refusal is a valid result. The scalar gate never trusts SVD error estimates.
The overlap premise is the declared balanced reflection-positive Hubbard
sector, not an arbitrary fermionic input.
"""
import argparse,json,sys,time
from fractions import Fraction as F
from pathlib import Path
from research.constructive_compression_20260916.model import fixture,mpo,exact_seed_cert
from research.constructive_compression_20260916.tensor_exact import state,inner,sqrt_upper,recurrence_residual
from research.constructive_compression_20260916.chebyshev_scalar import gate
from research.correlated_pair_20260913.mps_exact import State
from research.molecular_collective_20260913.core import digest


def check(proof):
    if proof['kind']!='hubbard_trace_chebyshev_proposal_v1':raise ValueError('Proof kind')
    start=time.monotonic();data=fixture(proof['rungs']);op=mpo(proof['rungs']);certs=proof['states']
    if len(certs)<2:raise ValueError('At least one recurrence step')
    if certs[0]!=exact_seed_cert(proof['rungs'],certs[0]['denominator']):raise ValueError('Exact trace seed required')
    b,ell,z=map(F,(proof['b'],proof['ell'],proof['z']));history=[];errors=[];older=None;previous=None;stats={}
    for i,c in enumerate(certs):
        # Use the existing strict schema/sector checker. Its norm is also exact.
        validated=State(data,c)
        if c['spin_counts']!=data['spin_counts']:raise ValueError('Balanced sector')
        v=state(c)
        n=F(validated.norm_integer,c['denominator']**(2*c['modes']))
        if i:
            rr,eta=recurrence_residual(v,previous,older,op,b,ell,stats);errors.append(eta)
            result=gate(b,ell,z,sqrt_upper(n),errors)
            error_floor=2*sum((z**(j+1)*e for j,e in enumerate(errors)),F(0))/(1-z*z)
            row={'step':i,'residual_squared_norm':str(rr),'residual_norm_upper':str(eta),
                 'residual_norm_upper_float':float(eta),'norm_upper':str(sqrt_upper(n)),
                 'score_float':float(result['score']),'error_floor_float':float(error_floor)}
            history.append(row);print(json.dumps({k:x for k,x in row.items() if not k.endswith('squared_norm')}),flush=True)
            if error_floor>=1:break
        older,previous=previous,v
    result={k:str(v) if isinstance(v,F) else v for k,v in result.items()}
    result.update({'status':'certified_lower' if result['accepted'] else 'criterion_not_met',
       'fixture_sha256':digest(data),'proof_sha256':digest(proof),'units':'t','history':history,'stats':stats,
       'declared_steps':len(certs)-1,'checked_steps':len(history),'enumerated_configurations':0,
       'global_matrix_entries':0,'replay_seconds':time.monotonic()-start,
       'all_declared_steps_replayed':len(history)==len(certs)-1,
       'scope':'balanced transformed 2xr Hubbard with U/t=8; PSD ground-overlap theorem'})
    if result['accepted'] and not result['all_declared_steps_replayed']:raise AssertionError('Incomplete accepted witness')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('proof',type=Path);a=p.parse_args()
    result=check(json.loads(a.proof.read_text()))
    if any(x in sys.modules for x in ('numpy','scipy','quimb','pyscf')):raise AssertionError('Numerical import in accepting path')
    a.proof.with_name(a.proof.stem+'_checked.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='history'},indent=2))
