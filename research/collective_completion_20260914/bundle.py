"""Independent two-sided acceptance of a new compact molecular proof."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
from research.collective_completion_20260914.spin_replay import check as balanced_check
from research.collective_completion_20260914.spin_screen import check as screened_check
from research.reconstruction_compression_20260914.upper_interval import check as upper_check

def run(fixture,state,certificate,upper_receipt,output,nonsinglet=None):
    started=time.monotonic();paths=[fixture,state,certificate,upper_receipt]+([nonsinglet] if nonsinglet else [])
    read=lambda p:json.loads(Path(p).read_text())
    data=read(fixture);psi=read(state);cert=read(certificate)
    proposed=F(read(upper_receipt)['upper_Ha']);grid=2**24;u=F(-((-proposed.numerator*grid)//proposed.denominator),grid)
    # The source receipt merely proposes u. Acceptance verifies the actual
    # rational MPS and uI-H anew, without trusting that receipt's bound.
    upper=upper_check(data,psi,u,True)
    if upper['status'] not in ('certified_upper','certified_upper_exact_fallback'):raise ValueError('Upper endpoint not certified')
    lower=screened_check(data,cert,read(nonsinglet)) if nonsinglet else balanced_check(data,cert)
    L=F(lower['lower'])
    if L>u:raise ValueError('Inconsistent endpoints')
    result={'lower_Ha':str(L),'upper_Ha':str(u),'width_Ha':str(u-L),'width_mHa':float(1000*(u-L)),
        'target_1p6mHa_met':u-L<=F(1,625),'upper':upper,'lower':lower,
        'source_endpoint_allowance_Ha':str(u-proposed),'total_replay_seconds':time.monotonic()-started,
        'all_input_sha256':{str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths},
        'lower_full_determinant_enumeration':False,'inherited_full_cubic_proof_used':False}
    with Path(output).open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ('upper','lower','all_input_sha256','upper_Ha','width_Ha')},indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('certificate');p.add_argument('upper_receipt');p.add_argument('output');p.add_argument('--nonsinglet');a=p.parse_args()
    run(a.fixture,a.state,a.certificate,a.upper_receipt,a.output,a.nonsinglet)
