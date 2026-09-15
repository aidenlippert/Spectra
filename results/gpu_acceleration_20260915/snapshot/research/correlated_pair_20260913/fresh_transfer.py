"""Freeze Hamiltonian-only builders before generating an unseen geometry."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

BASE=Path('results/correlated_pair_20260913')
SOURCES=['mps_spatial.py','mps_direct.py','mps_round.py','mps_exact.py','self_consistent/fixed_guide.py']


def hashes():
    return {p:hashlib.sha256((Path('research/correlated_pair_20260913')/p).read_bytes()).hexdigest() for p in SOURCES}


def freeze():
    rule={'created_unix_seconds':time.time(),'source_sha256':hashes(),'known_H6_spacings_A':[1.2,1.4,1.6,1.73,1.8],
          'new_H6_spacing_A':'191/100','upper':'Random spatial MPS seed20260913, initial bond4, sweep bonds8/16/48, ten RL sweeps, local tolerance1e-5, number penalty4, rational denominator1e9, exact charge projection and contraction.',
          'lower':'Fixed PH guide at chemical midpoint, five self-consistent iterations, damping1/2, coefficient denominator1e7, exact paired positive factors and L1 remainder.',
          'selection_inputs':'Hamiltonian and physical counts only. No FCI or cubic teacher; no reference accuracy used for stopping.',
          'required_output':'Full exact interval, even if wider than1.6mHa; no strong inherited lower exists on the new case.'}
    with (BASE/'transfer_freeze.json').open('x') as f:json.dump(rule,f,indent=2)
    print(json.dumps(rule,indent=2))


def generate():
    rule=json.loads((BASE/'transfer_freeze.json').read_text())
    if rule['source_sha256']!=hashes():raise ValueError('Frozen builder source changed')
    out=BASE/'fresh_h6_1p91'
    if out.exists():raise ValueError('Fresh fixture already exists')
    from research.mechanism_transfer_20260913.fixtures import build
    receipt=build(6,F(191,100),out,fci_control=False)
    receipt.update(FCI_reference_built=False,full_cubic_teacher_built=False,frozen_source_hashes_verified=True)
    (out/'generation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    if sys.argv[1:] == ['freeze']:freeze()
    elif sys.argv[1:] == ['generate']:generate()
    else:raise ValueError('freeze or generate required')
