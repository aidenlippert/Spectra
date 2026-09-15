"""Integrate a newly contracted MPS upper with an inherited cubic lower."""
from fractions import Fraction as F
import argparse
import json
from pathlib import Path
import sys
import time
from research.correlated_pair_20260913.mps_exact import check
from research.global_response_20260913.reference_diagnostic import lower
from research.certificate_scaling.streaming_reference_upper import upper


def run(case,folder):
    start=time.monotonic();data=json.loads(Path(f'results/certificate_scaling/active_space_ladder/{case}/fixture.json').read_text())
    cert=json.loads((folder/'state.json').read_text());uc=check(data,cert);U=F(uc['upper_Ha']);L,lc=lower(data,case)
    if U<L:raise ValueError('Inconsistent two-sided endpoints')
    # Explicit, separately identified comparison; removing this block leaves
    # the valid new interval unchanged. Its amplitudes never teach the MPS.
    meta=json.loads(Path(f'results/certificate_scaling/cubic_precision/intervals/final_{case}.json').read_text())
    ref=json.loads(Path(meta['reference']).read_text());R,rc=upper(data,ref['independent_upper'])
    out={'case':case,'upper':uc,'inherited_lower':lc,'lower_Ha':str(L),'upper_Ha':str(U),'width_mHa':float(1000*(U-L)),
         'width_Ha':str(U-L),'target_1p6mHa_met':U-L<=F(16,10000),
         'reference_upper_Ha':str(R),'upper_deterioration_mHa':float(1000*(U-R)),
         'upper_0p5mHa_target_met':U-R<=F(5,10000),'reference_comparison_replay':rc,
         'old_lower_discovery_still_inherited':True,'FCI_upper_dependency_for_new_interval':False,
         'reference_amplitudes_used_only_for_separate_comparison':True,'wall_seconds':time.monotonic()-start,
         'scope':'Exact supplied rational Hamiltonian in fixed N, with an M_S=0 variational upper; no physical-model uncertainty certified.'}
    if any(x in sys.modules for x in ('numpy','scipy','pyscf','quimb','cvxpy')):raise AssertionError('Numerical dependency during exact replay')
    (folder/'interval.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('upper','inherited_lower','reference_comparison_replay')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',choices=['h6','h8']);p.add_argument('folder',type=Path);a=p.parse_args();run(a.case,a.folder)
