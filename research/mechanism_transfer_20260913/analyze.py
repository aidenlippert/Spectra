"""Train on the frozen winning proof, then bind a rule before transfer."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import time
import numpy as np

from experiments.marginal_symbolic import canonical, mono
from research.trace_pricing_20260913.discovery import Model
from research.mechanism_transfer_20260913.rule import generate

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/mechanism_transfer_20260913'


def save(path, obj): path.write_text(json.dumps(obj,indent=2)+'\n')


def run():
    start = time.monotonic(); prior=ROOT/'results/trace_pricing_20260913'
    path=prior/'manifest.json'; old=json.loads(path.read_text())['files']
    bad=[n for n,v in old.items() if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=v['sha256']]
    if bad: raise ValueError('Frozen inputs changed: '+repr(bad))
    source=prior/'campaign/continue/round_10_batch_1/certificate.json'
    base=ROOT/'results/molecular_collective_20260913/campaign/h6'
    data=json.loads((base/'fixture.json').read_text()); tail=json.loads((base/'rank_10/tail.json').read_text())
    cert=json.loads(source.read_text()); model=Model(data,tail)
    blocks={b['name']:b for b in cert['anti_blocks']}; proof=[]; spectra=[]; monomial_factors=[]
    for gid,frame in enumerate(model.frames):
        b=blocks[frame['name']]
        D=np.array(b['directions'],float)/b['direction_denominator']
        R=np.array(b['factor'],float)/cert['denominator']
        W=model.frame_coefficients[gid]@D.T@R.T
        monomial_factors.append(W)
        chol=np.linalg.cholesky(model.monomial_traces[gid]).T
        Z=chol@W; proof.append((Z,chol)); values=np.linalg.svd(Z,compute_uv=False)**2
        cdf=np.cumsum(values)/sum(values)
        spectra.append({'group':gid,'factor_rows':len(R),'trace_size':float(sum(values)),
            'rank_for_99_percent_trace_size':int(np.searchsorted(cdf,.99)+1),
            'rank_for_999_percent_trace_size':int(np.searchsorted(cdf,.999)+1),
            'singular_values_squared':values.tolist()})
    flips=[]
    for gid,frame in enumerate(model.frames):
        words=model.spin_prepared[gid]['words']; mapped=[]
        for w in words:
            q=canonical(mono(tuple((c,i^1) for c,i in w)))
            if len(q)!=1: raise AssertionError('Spin flip must map a monomial to one monomial')
            mapped.append(next(iter(q.items())))
        target=next(j for j,f in enumerate(model.spin_prepared) if set(w for w,_ in mapped)==set(f['words']))
        if gid>=target: continue
        lookup={w:i for i,w in enumerate(model.spin_prepared[target]['words'])}
        P=np.zeros((len(mapped),len(mapped)))
        for j,(w,sign) in enumerate(mapped): P[lookup[w],j]=float(sign)
        a=P@monomial_factors[gid]; b=monomial_factors[target]
        ka=a@a.T; kb=b@b.T
        flips.append({'groups':[gid,target], 'relative_operator_Gram_difference':float(np.linalg.norm(ka-kb)/max(np.linalg.norm(ka),np.linalg.norm(kb)))})
    candidates=[]
    for powers in ((0.,),(.5,),(1.,),(0.,.5),(0.,1.),(.5,1.)):
        span,stats=generate(model,powers); total=0.; missed=0.; by_group=[]
        for gid,(Z,chol) in enumerate(proof):
            directions=[v['vector'] for v in span if v['group']==gid]
            C=model.frame_coefficients[gid]
            A=chol@C@(np.array(directions,float)/10**10).T if directions else np.empty((C.shape[0],0))
            Q=np.linalg.qr(A,mode='reduced')[0]; residual=Z-Q@(Q.T@Z)
            mass=float(np.sum(Z*Z)); error=float(np.sum(residual*residual))
            total+=mass; missed+=error
            by_group.append({'group':gid,'coverage':1-error/mass,'dimension':len(directions)})
        candidates.append({'powers':list(powers),'trace_factor_coverage':1-missed/total,'groups':by_group,'construction':stats})
    chosen=max((c for c in candidates if len(c['powers'])==2),key=lambda c:c['trace_factor_coverage'])
    entries=[abs(c) for b in cert['base_blocks']+cert['anti_blocks'] for row in b['factor'] for c in row if c]
    dirs=[abs(c) for b in cert['anti_blocks'] for row in b['directions'] for c in row if c]
    result={'source_certificate_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'frozen_manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'preserved_files':len(old),
        'candidates':candidates,'factor_spectra':spectra,'spin_flip_comparisons':flips,
        'factor_denominator':cert['denominator'],'direction_denominator':10**10,
        'maximum_factor_integer_bits':max(c.bit_length() for c in entries),
        'maximum_direction_integer_bits':max(c.bit_length() for c in dirs),
        'wall_seconds':time.monotonic()-start,
        'scope':'Training proof structure diagnostics only. Projection coverage is not retained energy-bound gain.'}
    save(OUT/'training_analysis.json',result)
    choice={'powers':chosen['powers'],'training_trace_factor_coverage':chosen['trace_factor_coverage'],
        'source_certificate_sha256':result['source_certificate_sha256'],
        'rule_source_sha256':hashlib.sha256((Path(__file__).parent/'rule.py').read_bytes()).hexdigest(),
        'training_analysis_sha256':hashlib.sha256((OUT/'training_analysis.json').read_bytes()).hexdigest(),
        'scope':'Frozen before held-out fixture generation or any new energy solve.'}
    save(OUT/'frozen_rule.json',choice); print(json.dumps(result),flush=True)


if __name__=='__main__': run()
