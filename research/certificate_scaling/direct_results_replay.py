"""Independent stdlib replay of direct-discovery molecular energy intervals.

Upper witnesses are supplied only after discovery has finished. Existing H4
witness lists enumerate 70 states; H6's FCI-derived control has 200 nonzero
states. Their generation is a validation cost, not scalable discovery.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import verify
from experiments.marginal_determinant_tree import DeterminantOracle

ROOT=Path(__file__).resolve().parents[2]

def replay_interval(cert):
    lower=F(verify(cert)['lower'])
    upper=DeterminantOracle(cert).upper(cert['sparse_independent_upper'])
    if upper<lower:raise ValueError('Inconsistent interval')
    return {'lower':str(lower),'upper':str(upper),'width':str(upper-lower),
            'lower_float':float(lower),'upper_float':float(upper),'width_float':float(upper-lower),
            'passes_0_0015_Ha':upper-lower<=F(3,2000),
            'upper_support':len(cert['sparse_independent_upper']['states'])}

def witness(name):
    if name=='h6':
        return json.loads((ROOT/'results/marginal_h6/reference_upper.json').read_text())['independent_upper']
    path=ROOT/('results/marginal_molecule_stress/accepted_degree3_certificate.json' if name=='square'
               else 'results/marginal_molecule/h4_compact_certificate.json')
    source=json.loads(path.read_text()); upper=source['independent_upper']
    states=[s for s in range(1<<source['modes']) if s.bit_count()==source['particles']]
    if len(states)!=len(upper['amplitudes']):raise ValueError('Upper witness order mismatch')
    return {'states':states,'amplitudes':upper['amplitudes']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True); rows=[]
    for path in sorted(a.input.rglob('certificate.json')):
        name=path.parent.name.split('_')[0]
        if name not in ('square','rectangle','h6'):continue
        cert=json.loads(path.read_text()); cert['sparse_independent_upper']=witness(name)
        start=time.monotonic(); receipt=replay_interval(cert)
        receipt['replay_seconds']=time.monotonic()-start
        label='_'.join(path.relative_to(a.input).parts[:-1])
        output=a.output/(label+'_interval.json')
        output.write_text(json.dumps(cert,separators=(',',':'))+'\n')
        receipt.update({'label':label,'certificate':str(output),'bytes':output.stat().st_size,
                        'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                        'certificate_sha256':hashlib.sha256(output.read_bytes()).hexdigest()})
        rows.append(receipt)
    (a.output/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps({'replayed':len(rows),'passing':sum(r['passes_0_0015_Ha'] for r in rows),
                      'seconds':sum(r['replay_seconds'] for r in rows)}))

if __name__=='__main__':main()
