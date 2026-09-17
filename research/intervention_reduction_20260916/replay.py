"""Standalone exact replay; refuse numerical modules in the accepting process."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import resource
import sys
from research.intervention_reduction_20260916.exact import check


def run(fixture, proposal, output, horizon=F(10)):
    if output.exists():
        raise FileExistsError(output)
    data=json.loads(fixture.read_text()); candidate=json.loads(proposal.read_text())
    receipt,moments=check(data,candidate,horizon=horizon)
    forbidden=[s for s in ('numpy','scipy','quimb','pyscf','numba') if s in sys.modules]
    if forbidden:
        raise AssertionError(('Numerical accepting dependency',forbidden))
    receipt['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    receipt['accepting_numerical_imports']=forbidden
    receipt['exact_reduced_model']={
        'metric':[[str(x) for x in row] for row in moments['metric']],
        'generators':[[[str(x) for x in row] for row in a] for a in moments['generators']]}
    output.write_text(json.dumps(receipt,separators=(',',':'))+'\n')
    print(json.dumps({k:receipt[k] for k in ('target_met','reduced_dimension','state_vector_error_float',
        'uniform_frobenius_error_float','selected_configurations','reached_configurations','initial_energy_float_Ha',
        'exact_replay_seconds','peak_RSS_bytes')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--horizon',default='10');a=p.parse_args()
    run(a.fixture,a.proposal,a.output,F(a.horizon))
