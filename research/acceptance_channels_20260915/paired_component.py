"""Time only the actual paired H12 component; this is not an energy replay."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import encode,expand_squares as original_expand
from research.acceptance_channels_20260915.paired_exact import expand_squares,CALLS,supported,validated


def run(certificate,output,original=False):
    started=time.monotonic();raw=certificate.read_bytes();core=json.loads(raw)['core']
    blocks=core['blocks'];selected=[];i=0
    while i<len(blocks)-1:
        words,factors=validated(blocks[i],core['denominator'],core['modes'],3)
        next_words,next_factors=validated(blocks[i+1],core['denominator'],core['modes'],3)
        if supported(words,core['modes']) and next_words==[tuple((1-c,k) for c,k in reversed(w)) for w in words] and next_factors==factors:
            selected.extend(blocks[i:i+2]);i+=2
        else:i+=1
    if not selected:raise ValueError('No strictly matching paired component')
    parse_seconds=time.monotonic()-started;begin=time.monotonic()
    evaluate=original_expand if original else expand_squares
    polynomial,stats=evaluate(selected,core['denominator'],core['modes'])
    elapsed=time.monotonic()-begin
    encoded=json.dumps(encode(polynomial),separators=(',',':')).encode()
    result={'kind':'exact_molecular_paired_component_only','source_certificate_sha256':hashlib.sha256(raw).hexdigest(),
        'backend':'original_CAR_expansion' if original else 'exact_quartic_contractions',
        'parse_and_select_seconds':parse_seconds,'evaluation_seconds':elapsed,'total_seconds':time.monotonic()-started,
        'selected_block_count':len(selected),'polynomial_terms':len(polynomial),'maximum_degree':max(map(len,polynomial),default=0),
        'polynomial_sha256':hashlib.sha256(encoded).hexdigest(),'stats':stats,'calls':CALLS,
        'complete_energy_certificate_replayed':False,'speedup_vs_original_measured':False}
    with output.open('x') as stream:json.dump(result,stream,indent=2)
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--original',action='store_true')
    a=p.parse_args();run(a.certificate,a.output,a.original)
