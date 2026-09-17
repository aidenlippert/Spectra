"""Require identical inputs and exact accepted endpoints for the verifier comparison."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def run(original,contracted,output):
    read=lambda p:json.loads(p.read_text())
    a,b=read(original/'original_interval.json'),read(contracted/'original_interval.json')
    for key in ('lower_Ha','upper_Ha','width_Ha'):
        if Fraction(a[key])!=Fraction(b[key]):raise ValueError(('Exact endpoints disagree',key))
    for name in ('fixture.json','nonsinglet.json','certificate.json'):
        if hashlib.sha256((original/'exact'/name).read_bytes()).digest()!=hashlib.sha256((contracted/'exact'/name).read_bytes()).digest():
            raise ValueError(('Verifier inputs differ',name))
    if hashlib.sha256((original/'mps/state.json').read_bytes()).digest()!=hashlib.sha256((contracted/'mps/state.json').read_bytes()).digest():
        raise ValueError('Upper witnesses differ')
    lower_a,lower_b=read(original/'exact/lower.json'),read(contracted/'exact/lower.json')
    if Fraction(lower_a['singlet']['residual_l1'])!=Fraction(lower_b['singlet']['residual_l1']):
        raise ValueError('Exact residual allowances disagree')
    for sector in ('singlet','all_nonsinglets'):
        for key in ('factor_rows','factor_nonzeros','square_polynomial_terms','residual_terms','residual_max_degree','mode_count','particle_number'):
            if lower_a[sector][key]!=lower_b[sector][key]:raise ValueError(('Proof statistics differ',sector,key))
    method=read(contracted/'paired_expansion.json')
    if method['forbidden_imports'] or not any(c['exact_matched_block_pairs'] for c in method['calls']):
        raise ValueError('The contracted standard-library path was not actually exercised')
    result={'original_replay':str(original),'contracted_replay':str(contracted),
        'identical_exact_endpoints':True,'identical_input_witnesses':True,'identical_exact_singlet_residual_allowance':True,
        'identical_factor_and_polynomial_statistics':True,
        'original_complete_seconds':a['complete_seconds'],'contracted_complete_seconds':b['complete_seconds'],
        'original_time_divided_by_contracted_time':a['complete_seconds']/b['complete_seconds'],
        'width_mHa':a['width_mHa'],'paired_calls':method['calls'],
        'comparison_scope':'One complete replay each on this host. Excludes discovery and unsuccessful earlier replay; not an end-to-end solver speedup.'}
    with output.open('x') as stream:json.dump(result,stream,indent=2)
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('original',type=Path);p.add_argument('contracted',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();run(a.original,a.contracted,a.output)
