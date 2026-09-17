"""Fresh exact kernel construction followed by reusable small-matrix queries."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import resource
import time
from research.intervention_reduction_20260916.kernel import CheckedKernel
from research.intervention_reduction_20260916.trajectory import check_trajectory
from research.intervention_reduction_20260916.robust import extend_control_neighborhood


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('output',type=Path);p.add_argument('trajectories',type=Path,nargs='+')
    p.add_argument('--max-deviation-Ha',default='1/400000');a=p.parse_args()
    a.output.mkdir(exist_ok=False,parents=True);started=time.monotonic()
    data=json.loads(a.fixture.read_text());proposal=json.loads(a.proposal.read_text())
    trajectories=[json.loads(path.read_text()) for path in a.trajectories]
    kernel=CheckedKernel(data,proposal,trajectories[0]['observed_spatial_orbital'])
    del data,proposal
    queries=[]
    for i,(path,trajectory) in enumerate(zip(a.trajectories,trajectories)):
        query_start=time.monotonic();nominal=check_trajectory(None,None,trajectory,_kernel=kernel)
        budgets=[F(a.max_deviation_Ha)*F(trajectory['horizon_atomic_time'])]*len(kernel.query_proposal['controls'])
        result=extend_control_neighborhood(kernel.query_data,kernel.query_proposal,trajectory,nominal,budgets)
        result['compiled_query_seconds']=time.monotonic()-query_start
        result['kernel_rebuilt_once_from_original_inputs']=True
        (a.output/f'query_{i+1}.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
        queries.append({'input':str(path),'output':f'query_{i+1}.json','seconds':result['compiled_query_seconds'],
            'target_met':result['target_met'],'uniform_in_time_target_met':result['uniform_in_time_target_met'],
            'state_bound':result['normalized_state_error_float'],
            'population_change':result['population_change_interval_float']})
    report={'kernel':kernel.receipt,'queries':queries,'total_seconds':time.monotonic()-started,
        'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'claim':'query matrices are small after charged exact construction; no enumeration-free discovery claim'}
    (a.output/'bundle.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
