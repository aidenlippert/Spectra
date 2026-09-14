"""New accepting caller: rigorous fast upper plus preserved strong lower."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from research.reconstruction_compression_20260914.inputs import OUT,dump
from research.reconstruction_compression_20260914.upper_interval import check
from research.global_response_20260913.reference_diagnostic import lower

def run(case):
    start=time.monotonic();inputs=json.loads((OUT/'frozen_inputs.json').read_text())[case]
    data=json.loads(Path(inputs['fixture']).read_text());state=json.loads(Path(inputs['state']).read_text())
    comparison=json.loads((OUT/'upper'/case/'comparison.json').read_text());u=F(comparison['shared_comparison_endpoint_Ha'])
    upper=check(data,state,u,True)
    if upper['status'] not in ('certified_upper','certified_upper_exact_fallback'):raise ValueError('Upper endpoint not certified')
    L,proof=lower(data,case)
    if L!=F(inputs['lower_Ha']) or L>u:raise ValueError('Lower endpoint mismatch')
    rec={'case':case,'lower_Ha':str(L),'upper_Ha':str(u),'width_mHa':float(1000*(u-L)),
         'target_1p6mHa_met':u-L<=F(1,625),'upper':upper,'lower':proof,'seconds':time.monotonic()-start,
         'compressed_lower_achieved':False,'inherited_full_cubic_discovery':True,
         'original_exact_upper_Ha':inputs['upper_Ha'],'charged_upper_allowance_Ha':str(u-F(inputs['upper_Ha']))}
    folder=OUT/'bundles';folder.mkdir(exist_ok=True);dump(folder/f'{case}.json',rec)
    print(json.dumps({k:v for k,v in rec.items() if k not in ('upper','lower','original_exact_upper_Ha','charged_upper_allowance_Ha')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',choices=['h6','h8']);run(p.parse_args().case)
