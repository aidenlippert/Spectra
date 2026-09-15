"""Run independent library probes sequentially and retain every process outcome."""
import json
from pathlib import Path
import subprocess
import sys
import time


def run():
    root=Path(__file__).resolve().parents[2]
    out=root/'results/nvidia_followup_20260915'
    if (out/'probe_jobs.json').exists():
        raise ValueError('Preserve existing probe jobs')
    prefix=[sys.executable,'-B','-m']
    cases=root/'results/transfer_solver_20260915'
    bench='research.nvidia_followup_20260915.library_bench'
    exact='research.nvidia_followup_20260915.flint_probe'
    jobs=[('sparse_quotient_tests',60,prefix+['unittest','research.nvidia_followup_20260915.test_sparse_quotient']),
        ('flint_gram_tests',60,prefix+['unittest','research.nvidia_followup_20260915.test_flint_squares']),
        ('h8_sparse',180,prefix+[bench,'sparse',str(cases/'cases/h8_cold'),str(out/'h8_sparse.json')]),
        ('h10_sparse',180,prefix+[bench,'sparse',str(cases/'adaptive/h10_correlated_guide'),str(out/'h10_sparse.json')]),
        ('h10_tensor',180,prefix+[bench,'tensor',str(cases/'rotated_h10'),str(out/'h10_tensor.json')]),
        ('h8_stdlib',240,prefix+[exact,str(cases/'cases/h8_cold'),str(out/'h8_stdlib.json'),'--stdlib']),
        ('h8_flint',240,prefix+[exact,str(cases/'cases/h8_cold'),str(out/'h8_flint.json')]),
        ('h8_flint_gram',240,prefix+[exact,str(cases/'cases/h8_cold'),str(out/'h8_flint_gram.json'),'--grams'])]
    results=[]
    for name,seconds,command in jobs:
        start=time.monotonic()
        with (out/(name+'.log')).open('x') as stream:
            p=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'_time.txt')),
                'timeout','--signal=TERM','--kill-after=10s',str(seconds)]+command,
                cwd=root,stdout=stream,stderr=subprocess.STDOUT)
        row={'name':name,'seconds':time.monotonic()-start,'exit_code':p.returncode,
            'status':'passed' if p.returncode==0 else 'failed','command':command}
        results.append(row)
        (out/'probe_jobs.json').write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(row),flush=True)


if __name__=='__main__':
    run()
