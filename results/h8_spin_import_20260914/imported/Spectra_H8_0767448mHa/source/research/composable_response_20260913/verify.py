"""Read-only exact replay, including the preserved numerical milestone."""
import argparse
import hashlib
import json
import sys
import time

from research.compact_response_20260913 import program,closure
from research.composable_response_20260913 import joint,recursive,finite_controls
from research.ch2_validation_20260913 import certify

ROOT=program.ROOT
PARENT_SHA='619ee08d61d7260373d1b364d87cfe533c7ba94bc36f0658bc3a98431682cde3'


def integrity(path,expected_hash=None):
    raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    if expected_hash is not None and sha!=expected_hash:raise ValueError('Changed parent manifest')
    entries=json.loads(raw)['files']
    for name,expected in entries.items():
        data=(ROOT/name).read_bytes()
        if len(data)!=expected['bytes'] or hashlib.sha256(data).hexdigest()!=expected['sha256']:
            raise ValueError('Manifest mismatch: '+name)
    return {'files_verified':len(entries),'manifest_sha256':sha}


def replay(allow_unsealed=False):
    start=time.monotonic();parent=integrity(program.OUT/'manifest.json',PARENT_SHA)
    path=joint.OUT/'manifest.json'
    current=integrity(path) if path.exists() else None
    if current is None and not allow_unsealed:raise ValueError('The new result has not been sealed')
    data,tail,reference=program.load_case('h6')
    first=json.loads((program.OUT/'h6_program.json').read_text())
    milestone=closure.check(data,tail,first,json.loads((program.OUT/'expanded_closure_certificate.json').read_text()),reference)
    finite=finite_controls.replay();composed=[];obstructions=[]
    for name in ('h6','fresh_h6_1p6','h8'):
        data,tail,_,first=joint.load_case(name)
        composed.append({'case':name,'receipt':joint.check(data,tail,first,json.loads((joint.OUT/f'{name}_joint_response.json').read_text()))})
        if name!='fresh_h6_1p6':
            r=recursive.check(data,tail,first,json.loads((joint.OUT/f'{name}_recursive.json').read_text()))
            obstructions.append({'case':name,'receipt':r})
    data=json.loads((certify.OUT/'fixture.json').read_text())
    spin=[certify.check(data,json.loads((certify.OUT/f'spin_{s}_certificate.json').read_text())) for s in (0,1)]
    low=certify.F(spin[1]['lower_Ha'])-certify.F(spin[0]['upper_Ha'])
    high=certify.F(spin[1]['upper_Ha'])-certify.F(spin[0]['lower_Ha'])
    forbidden=[x for x in ('numpy','scipy','cvxpy','pyscf') if x in sys.modules]
    if forbidden:raise AssertionError('Numerical import during exact replay')
    return {'parent_integrity':parent,'current_integrity':current,'preserved_H6':milestone,
        'matched_controls_and_transfer':finite,'composed_responses':composed,'scalar_obstructions':obstructions,
        'CH2':{'states':spin,'gap_convention':'E_T-E_S','gap_lower_Ha':str(low),'gap_upper_Ha':str(high),
            'gap_width_mHa':float(1000*(high-low))},
        'numerical_packages_loaded':forbidden,'wall_seconds':time.monotonic()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--allow-unsealed',action='store_true');args=parser.parse_args()
    print(json.dumps(replay(args.allow_unsealed),indent=2))
