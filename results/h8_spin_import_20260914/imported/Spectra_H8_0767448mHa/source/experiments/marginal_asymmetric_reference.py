"""Controlled full-quartic reference with an external process deadline."""
from fractions import Fraction as F
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def model(epsilon,modes=10):
    from experiments.marginal_collective import hopping_polynomial
    from experiments.marginal_symmetry_transfer import hopping_perturbation
    from experiments.marginal_symbolic import add
    return add(hopping_polynomial(modes,F(1,5)),hopping_perturbation(modes,epsilon)[0])


def worker(epsilon,out):
    from experiments.marginal_dual_exact import full_word_blocks
    from experiments.marginal_asymmetric_adapt import exchange_subspaces
    from experiments.marginal_adaptive import assemble_and_solve
    from experiments.marginal_coefficient import export
    from experiments.marginal_symmetry_transfer import integer_upper
    from experiments.marginal_transfer_verify import replay
    modes=10;particles=5;h=model(epsilon);blocks=[];original=full_word_blocks(modes,4)
    for block in original:
        for j,tr in enumerate(exchange_subspaces(block['words'],modes)):
            blocks.append({'name':block['name']+f':lr{j}','words':block['words'],'basis_transform':tr})
    dims=[b['basis_transform'].shape[1] for b in blocks]
    structure={'modes':modes,'particles':particles,'epsilon':str(epsilon),'original_blocks':len(original),
               'original_maximum':max(len(b['words']) for b in original),'exchange_blocks':len(blocks),
               'exchange_maximum':max(dims),'psd_scalars':sum(d*(d+1)//2 for d in dims),
               'groups':[[i,i+5] for i in range(5)]}
    (out/'structure.json').write_text(json.dumps(structure,indent=2)+'\n');print(json.dumps(structure),flush=True)
    proposal,_=assemble_and_solve(h,modes,particles,blocks,degree=4,groups=structure['groups'],eps=1e-8,residual_penalty=True)
    numerical={k:proposal[k] for k in ('b','numerical_objective','status','build_seconds','solve_seconds')}
    (out/'numerical_proposal.json').write_text(json.dumps(numerical,indent=2)+'\n');print(json.dumps(numerical),flush=True)
    certificate,_=export(h,modes,particles,blocks,proposal,denominator=10**10,operator_degree=4)
    certificate['independent_upper']=integer_upper(h,modes,particles)
    (out/'certificate.json').write_text(json.dumps(certificate)+'\n')
    receipt=replay(certificate);receipt.update(numerical)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


def run(epsilon=F(1,1000),budget=120):
    if not 0<budget<=600:raise ValueError('Budget must be in (0,600] seconds')
    out=ROOT/'results/marginal_asymmetric_reference'/('controlled_'+str(epsilon).replace('/','_'))
    out.mkdir(parents=True,exist_ok=True)
    if (out/'process.json').exists():raise ValueError('Preserve previous controlled run; choose a new epsilon or archive it explicitly')
    started=time.monotonic();observations=[]
    with (out/'worker.log').open('w') as log:
        process=subprocess.Popen([sys.executable,'-m','experiments.marginal_asymmetric_reference','--worker','--epsilon',str(epsilon)],
                                 cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        def record(status):
            result={'pid':process.pid,'elapsed_seconds':time.monotonic()-started,'budget_seconds':budget,'status':status,'returncode':process.poll()}
            observations.append(result);(out/'process.json').write_text(json.dumps(observations,indent=2)+'\n');return result
        record('running')
        while process.poll() is None:
            remaining=budget-(time.monotonic()-started)
            if remaining<=0:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait(timeout=5)
                result=record('wall_budget_terminated');break
            try:process.wait(timeout=min(10,remaining))
            except subprocess.TimeoutExpired:record('running')
        else:result=record('completed' if process.returncode==0 else 'worker_failed')
    print(json.dumps(result),flush=True);return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--epsilon',default='1/1000');parser.add_argument('--budget',type=float,default=120);parser.add_argument('--worker',action='store_true')
    args=parser.parse_args();epsilon=F(args.epsilon)
    if args.worker:worker(epsilon,ROOT/'results/marginal_asymmetric_reference'/('controlled_'+str(epsilon).replace('/','_')))
    else:run(epsilon,args.budget)
