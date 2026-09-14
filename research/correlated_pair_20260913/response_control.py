"""Matched complete-bound control: does response help an inherited SOS path?"""
from fractions import Fraction as F
import json
from pathlib import Path
import time
from research.global_response_20260913 import global_program as gp
from research.global_response_20260913.reference_diagnostic import lower,scalar_terminal_margin


def run():
    rows=[]
    for case,folder in (('h6','h6_b48_real'),('h8','h8_spatial_warm144')):
        start=time.monotonic();data,tail,_=gp.load_case(case);L,lc=lower(data,case)
        # A common, slightly weakened lower gives identical width to both paths.
        b=L-F(2,10**6);build=time.monotonic();prior=json.loads((gp.OUT/f'{case}_global.json').read_text())
        bounds=gp.bounds(data,tail,b,prior['first_sector'],prior['joint_sector'])
        cert={**prior,'target_Ha':str(b),'budget_Ha':str(F(1,10**6)),'partition':bounds['partition'],
              'response':gp.scalar_program(bounds,F(1,10**6))}
        discovery={'sector_gap_proofs':'Inherited from accepted global_response campaign; exact replay retained.',
                   'new_work':'Scalar response construction at the matched target; no numerical proposer needed.'}
        rec=gp.check(data,tail,cert)
        if rec['status']!='certified_response':raise ValueError('Response unavailable')
        terminal=scalar_terminal_margin(L,b,[F(rec['residual_penalty_Ha'])])
        if not terminal['positive']:raise ValueError('Missing terminal positivity')
        extra=time.monotonic()-build
        interval=json.loads(Path(f'results/correlated_pair_20260913/mps/{folder}/interval.json').read_text());U=F(interval['upper_Ha'])
        rows.append({'case':case,'matched_lower_Ha':str(b),'matched_width_mHa':float(1000*(U-b)),
            'direct_path':'Verified inherited cubic lower L implies weaker b directly.',
            'response_path':'Same inherited global SOS implies Schur coercivity; response residual consumes part of its margin.',
            'terminal':terminal,'response':cert,'response_replay':rec,'response_discovery':discovery,
            'additional_response_discovery_and_replay_seconds':extra,'common_inherited_lower_replay_seconds':lc['replay_seconds'],
            'full_new_structured_terminal_proof':False,'result':'No accuracy or discovery gain in this matched inherited-SOS control.',
            'wall_seconds':time.monotonic()-start})
    out=Path('results/correlated_pair_20260913/response_comparison.json');out.write_text(json.dumps({'cases':rows,
        'limitation':'No competitive response-assisted structured lower was obtained. Operator action savings from the prior campaign are preserved but are not a solver speedup here.'},indent=2)+'\n');print(json.dumps([{'case':r['case'],'width_mHa':r['matched_width_mHa'],'extra_seconds':r['additional_response_discovery_and_replay_seconds']} for r in rows],indent=2))


if __name__=='__main__':run()
