"""Same rational state and declared endpoint; integer oracle cross-check."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from research.reconstruction_compression_20260914.inputs import OUT,dump
from research.correlated_pair_20260913.mps_exact import check as exact_check
from research.reconstruction_compression_20260914.upper_interval import check

def run(case):
    frozen=json.loads((OUT/'frozen_inputs.json').read_text())[case]
    data=json.loads(Path(frozen['fixture']).read_text());state=json.loads(Path(frozen['state']).read_text());U=F(frozen['upper_Ha'])
    # Declare one shared comparison endpoint before either checker runs. The
    # original exact U remains the endpoint used by lower-search comparisons.
    quantum=2**24;u=F(-((-U.numerator*quantum)//U.denominator),quantum)
    folder=OUT/'upper'/case;folder.mkdir(parents=True,exist_ok=False)
    approx=check(data,state,u);dump(folder/'enclosure.json',approx)
    exact=exact_check(data,state);dump(folder/'exact.json',exact)
    if F(exact['upper_Ha'])!=U:raise ValueError('Frozen MPS upper changed')
    norm=F(exact['norm_squared']);gap=(u-U)*norm
    if not F(approx['norm_enclosure'][0])<=norm<=F(approx['norm_enclosure'][1]):raise AssertionError('Norm enclosure missed exact value')
    if not F(approx['endpoint_inequality_enclosure'][0])<=gap<=F(approx['endpoint_inequality_enclosure'][1]):raise AssertionError('Energy enclosure missed exact value')
    at_equality=check(data,state,U);dump(folder/'exact_endpoint_attempt.json',at_equality)
    invalid=check(data,state,U-F(1,10000));dump(folder/'invalid_endpoint.json',invalid)
    if invalid['status'] not in ('endpoint_refuted','ambiguous_refused'):raise AssertionError('Invalid endpoint accepted')
    rec={'case':case,'same_state':frozen['state'],'exact_frozen_upper_Ha':str(U),'shared_comparison_endpoint_Ha':str(u),
         'charged_endpoint_allowance_Ha':str(u-U),'charged_endpoint_allowance_mHa':float(1000*(u-U)),
         'interval_checker_status':approx['status'],'exact_checker_accepts_same_endpoint':U<=u,
         'exact_endpoint_interval_status':at_equality['status'],'invalid_endpoint_status':invalid['status'],
         'exact_oracle_contains_all_enclosures':True,'integer_replay_seconds':exact['replay_seconds'],
         'interval_replay_seconds':approx['total_seconds'],'speed_ratio':exact['replay_seconds']/approx['total_seconds'],
         'original_interval_unchanged':True,'comparison_interval_mHa':float(1000*(u-F(frozen['lower_Ha'])))}
    dump(folder/'comparison.json',rec);print(json.dumps(rec,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',choices=['h6','h8']);run(p.parse_args().case)
