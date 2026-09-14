"""Freeze rules before generating the new geometry; no reference solve."""
from fractions import Fraction as F
import argparse
import json
from pathlib import Path
import subprocess
import time
from research.reconstruction_compression_20260914.inputs import ROOT,OUT,sha,dump

SOURCES=['research/reconstruction_compression_20260914/'+s for s in ('reduced.py','moments.py','ideal_rank.py','replay.py','transfer.py')]+[
    'research/correlated_pair_20260913/'+s for s in ('mps_spatial.py','mps_direct.py','mps_round.py','mps_exact.py')]+[
    'research/mechanism_transfer_20260913/fixtures.py','experiments/marginal_symbolic.py','experiments/marginal_coefficient.py',
    'research/certificate_scaling/adaptive_block_discovery.py','research/certificate_scaling/direct_sparse_discovery.py']

def hashes():return {p:sha(ROOT/p) for p in SOURCES}

def freeze():
    files=subprocess.run(['rg','--files','results','-g','*fixture.json'],check=True,capture_output=True,text=True).stdout.splitlines()
    known=[]
    for name in files:
        data=json.loads(Path(name).read_text())
        if data.get('natoms')==6 and data.get('kind')=='straight_hydrogen_chain_fixture_v1':
            z=data['geometry'][1][1][2]-data['geometry'][0][1][2]
            known.append({'path':name,'spacing_A':str(F(str(z)))})
    spacing=F(207,100)
    if any(F(k['spacing_A'])==spacing for k in known):raise ValueError('Chosen geometry already used')
    rule={'created_unix_seconds':time.time(),'source_sha256':hashes(),'known_chain_fixtures':known,
          'fresh_H6_spacing_A':str(spacing),'fixture_family':'Straight STO-3G RHF-canonical H6, rational coefficient grid1e12; FCI disabled',
          'MPS':'Independent spatial MPS, seed20260913, bonds8/16/48, ten RL sweeps, exact charge projection, rational denominator1e9',
          'lower_mode':'low','rank_cap':32,'solver_seconds':65,'exact_quartic_ideal_elimination':True,'solver':'SCS direct, epsilon1e-8, max30000iterations',
          'acceptance':'Existing full exact CAR checker and coefficient L1 remainder; full interval even if target fails',
          'teacher_or_reference_state_allowed':False,'upper_and_lower_must_be_constructed_on_fresh_input':True}
    with (OUT/'transfer_freeze.json').open('x') as f:json.dump(rule,f,indent=2)
    print(json.dumps({k:v for k,v in rule.items() if k not in ('source_sha256','known_chain_fixtures')},indent=2))

def generate():
    rule=json.loads((OUT/'transfer_freeze.json').read_text())
    if hashes()!=rule['source_sha256']:raise ValueError('Frozen source changed')
    folder=OUT/'fresh_h6_2p07'
    if folder.exists():raise ValueError('Fresh fixture already exists')
    from research.mechanism_transfer_20260913.fixtures import build
    result=build(6,F(rule['fresh_H6_spacing_A']),folder,fci_control=False)
    result.update(frozen_sources_verified=True,full_cubic_teacher_built=False,FCI_built=False)
    dump(folder/'generation_receipt.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['freeze','generate']);a=p.parse_args()
    freeze() if a.operation=='freeze' else generate()
