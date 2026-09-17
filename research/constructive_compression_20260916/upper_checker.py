"""Exact standard-library replay, with energy units explicitly t, not hartree."""
import argparse,json,sys
from pathlib import Path
from research.constructive_compression_20260916.model import fixture
from research.correlated_pair_20260913.mps_exact import check

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args()
    data=json.loads((a.directory/'fixture.json').read_text());cert=json.loads((a.directory/'state.json').read_text())
    if data!=fixture(data['rungs']):raise ValueError('Model binding')
    if cert['spin_counts']!=data['spin_counts']:raise ValueError('Balanced sector binding')
    result=check(data,cert)
    result['upper_t']=result.pop('upper_Ha');result['upper_float_t']=result.pop('upper_float_Ha')
    result['scope']='balanced spin projection of the transformed Hubbard model; no lower bound'
    if any(x in sys.modules for x in ('numpy','scipy','quimb','pyscf')):raise AssertionError('Numeric dependency')
    (a.directory/'upper_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
