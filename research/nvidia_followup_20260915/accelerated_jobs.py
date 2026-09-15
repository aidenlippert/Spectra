"""Iteration-matched comparisons, then a fresh integrated H10 proposal/replay."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def run():
    root=Path(__file__).resolve().parents[2]
    out=root/'results/nvidia_followup_20260915'
    if (out/'accelerated_jobs.json').exists():
        raise ValueError('Preserve existing accelerated runs')
    tests=json.loads((out/'probe_jobs.json').read_text())
    for name in ('sparse_quotient_tests','flint_gram_tests'):
        if not any(x['name']==name and x['status']=='passed' for x in tests):
            raise ValueError('Run focused arithmetic tests before integrated experiments')
    cases=root/'results/transfer_solver_20260915'
    h8=cases/'cases/h8_cold'
    h10=cases/'adaptive/h10_correlated_guide'
    prefix=[sys.executable,'-B','-m']
    solver='research.nvidia_followup_20260915.solve'
    base='research.gpu_acceleration_20260915.solve'
    bench='research.nvidia_followup_20260915.library_bench'
    jobs=[('h8_sparse_context_corrected',180,prefix+[bench,'sparse',str(h8),str(out/'h8_sparse_corrected.json')]),
        ('h10_sparse_context_corrected',180,prefix+[bench,'sparse',str(h10),str(out/'h10_sparse_corrected.json')]),
        ('h10_tensor_retry',180,prefix+[bench,'tensor',str(cases/'rotated_h10'),str(out/'h10_tensor_retry.json')])]
    for label,case in [('h8',h8),('h10',h10)]:
        jobs.append((label+'_dense_200',400,prefix+[base,str(case),'dense_200','--seconds','360',
            '--mu','2','--backend','cpu_evd','--iterations','200']))
        for backend,normal,tag in [('cpu_evd','superlu','sparse_cpu_200'),
            ('hybrid','superlu','sparse_hybrid_200'),('hybrid','cudss','sparse_cudss_200')]:
            jobs.append((label+'_'+tag,240,prefix+[solver,str(case),tag,'--seconds','210',
                '--backend',backend,'--normal',normal,'--iterations','200']))
    jobs.append(('h10_sparse_cudss_adaptive',900,prefix+[solver,str(h10),'sparse_cudss_adaptive',
        '--seconds','850','--backend','hybrid','--normal','cudss']))
    jobs.append(('h10_compiled_complete_replay',600,[sys.executable,'-B','-m','research.nvidia_followup_20260915.replay',
        str(h10),'sparse_cudss_adaptive',str(h10/'sparse_cudss_adaptive/exact_compiled'),
        '--rotated',str(cases/'rotated_h10'),'--compiled-rationals']))
    source_hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (root/'research/nvidia_followup_20260915').glob('*.py')}
    (out/'accelerated_protocol.json').write_text(json.dumps({'source_sha256':source_hashes,'jobs':jobs,
        'iteration_controls':200,'fixed_controls_mu':2.,'full_H10_start':'zero Gram and ideal coordinates, this pass correlated moments',
        'fresh_problem_preparation_costs_additional':True,'one_heavy_remote_process_at_a_time':True},indent=2)+'\n')
    results=[]
    for name,seconds,command in jobs:
        if name=='h10_compiled_complete_replay' and not (h10/'sparse_cudss_adaptive/export/certificate.json').exists():
            results.append({'name':name,'status':'blocked_by_missing_proposal','seconds':0.})
            continue
        start=time.monotonic()
        with (out/(name+'.log')).open('x') as stream:
            p=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'_time.txt')),
                'timeout','--signal=TERM','--kill-after=10s',str(seconds)]+command,
                cwd=root,stdout=stream,stderr=subprocess.STDOUT)
        row={'name':name,'seconds':time.monotonic()-start,'exit_code':p.returncode,
            'status':'passed' if p.returncode==0 else 'failed','command':command}
        results.append(row)
        (out/'accelerated_jobs.json').write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(row),flush=True)


if __name__=='__main__':
    run()
