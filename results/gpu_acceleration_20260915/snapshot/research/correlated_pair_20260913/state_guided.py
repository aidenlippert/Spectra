"""Let exact correlated moments propose factor corrections; retain by energy."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import decode,encode,product,add,scale,adj,mono
from research.correlated_pair_20260913.mps_exact import State
from research.correlated_pair_20260913.self_consistent.fixed_guide import ph,check,expand,rounded
from research.joint_patterns_20260913.core import anticommutator
from research.molecular_collective_20260913.core import digest


def run(data,statecert,basecert,out):
    start=time.monotonic();state=State(data,statecert);m,n=data['modes'],data['particles'];h=decode(data['hamiltonian'],m,4);U=state.expectation(h)
    weights=list(map(F,basecert['weights']));taus=[decode(t,m,3) for t in basecert['taus']];optimal=[]
    for i,t in enumerate(taus):
        c=add(product(mono(((1,i),)),t),product(adj(t),mono(((0,i),))))
        anti=anticommutator(t,t);a=state.expectation(ph(anti,n));b=state.expectation(ph(c,n))
        if a<0 or (a==0 and b!=0):raise ValueError('Inconsistent exact positive-factor moments')
        optimal.append(-b/(2*a) if a else F(1))
    rows=[];best=None;bestreceipt=None
    for blend in (F(0),F(1,4),F(1,2),F(1)):
        trial=[rounded(scale(t,1+blend*(x-1)),10**7) for t,x in zip(taus,optimal)]
        cert={**basecert,'taus':list(map(encode,trial))};rec=check(data,cert)
        sos,_=expand(trial,weights);G=state.expectation(ph(sos,n));eps=F(rec['residual_norm_bound_Ha']);width=U-F(rec['lower_Ha'])
        if not F(0)<=G<=width<=G+2*eps:raise AssertionError('State-proof energy identity inequality failed')
        rows.append({'blend':str(blend),'lower_Ha':rec['lower_Ha'],'width_mHa':float(1000*width),
                     'positive_factor_expectation_Ha':str(G),'residual_norm_bound_Ha':str(eps),'replay':rec})
        if bestreceipt is None or F(rec['lower_Ha'])>F(bestreceipt['lower_Ha']):best,bestreceipt=cert,rec
    result={'state_sha256':digest(statecert),'guide_certificate_sha256':digest(basecert),'optimal_factor_scales':list(map(str,optimal)),
            'scale_float':list(map(float,optimal)),'candidates':rows,'best_lower_Ha':bestreceipt['lower_Ha'],
            'selection':'Largest exactly verified lower, with the same exact MPS upper; factor expectation alone is not acceptance.',
            'FCI_or_full_cubic_teacher_used':False,'seconds':time.monotonic()-start}
    out.mkdir(parents=True,exist_ok=True);(out/'certificate.json').write_text(json.dumps(best,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'seconds':result['seconds'],'widths_mHa':[r['width_mHa'] for r in rows],'scales':result['scale_float']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('lower');p.add_argument('output',type=Path);a=p.parse_args()
    run(json.loads(Path(a.fixture).read_text()),json.loads(Path(a.state).read_text()),json.loads(Path(a.lower).read_text()),a.output)
