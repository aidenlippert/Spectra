"""Construction, independent acceptance, and adaptive refinement API."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from math import comb
import json,resource,subprocess,sys
from . import exact,io
from .solver import cocg
from .variational import sweep_solve,enrich

def solve(request,out,*,bonds=(48,96,160),seed_steps=12,sweeps=4,backend='gmp',verbose=True):
    """Solve a declared query; preserve every failure and accept only exact bounds."""
    start=perf_counter();out=Path(out)
    if out.exists():raise FileExistsError(out)
    if request.get('kind')!='tensor_hubbard_response_request_v1':raise ValueError('request kind')
    spec=request['model'];n,_,_=exact.validate_model(spec)
    if request.get('source')!=exact.source(spec):raise ValueError('source')
    frequencies=request.get('frequencies',[])
    if not 1<=len(frequencies)<=32:raise ValueError('frequency count')
    for row in frequencies:
        exact.rat(row['omega'])
        if exact.rat(row['eta'])<=0:raise ValueError('positive broadening')
    tolerance=exact.rat(request.get('target_radius','1/1000'))
    if tolerance<=0:raise ValueError('positive target')
    if not bonds or any(type(b) is not int or not 1<=b<=256 for b in bonds) or list(bonds)!=sorted(set(bonds)):raise ValueError('increasing bond limits up to 256')
    if type(seed_steps) is not int or seed_steps<1 or type(sweeps) is not int or sweeps<1:raise ValueError('iteration budget')
    if backend not in ('gmp','python'):raise ValueError('backend')
    out.mkdir(parents=True);io.dump(out/'request.json',request);results=[];warm=None
    for i,frequency in enumerate(frequencies):
        folder=out/f'query_{i:02d}';folder.mkdir();one_request=dict(request,frequencies=[frequency]);io.dump(folder/'request.json',one_request)
        z=complex(float(F(frequency['omega'])),float(F(frequency['eta'])));x=warm.copy() if warm is not None else None;phases=[];best=None
        for cap in bonds:
            if x is not None and x.bond>cap:continue
            began=perf_counter();phase=dict(bond_limit=cap)
            if x is None:x,phase['seed']=cocg(spec,z,max_bond=cap,maxiter=seed_steps,assess_every=seed_steps,rtol=.5*(float(tolerance)*z.imag)**.5,verbose=verbose)
            elif x.bond<cap:x,phase['growth']=enrich(spec,z,x,cap)
            x,phase['optimization']=sweep_solve(spec,z,x,sweeps=sweeps,target=.55*float(tolerance),verbose=verbose)
            candidate=io.candidate(spec,[z],[x],bits=30)
            for key in ('omega','eta'):candidate['queries'][0][key]=str(F(frequency[key]))
            name=f'bond_{cap:03d}';cp=folder/f'{name}.candidate.json';rp=folder/f'{name}.receipt.json';io.dump(cp,candidate)
            verify_start=perf_counter();cmd=[sys.executable,'-B','-S','-m','spectra_tensor.exact',str(cp.resolve()),'--request',str((folder/'request.json').resolve()),'--out',str(rp.resolve()),'--backend',backend]
            process=subprocess.run(cmd,capture_output=True,text=True,cwd=Path(__file__).resolve().parent.parent)
            (folder/f'{name}.checker.log').write_text(process.stdout+process.stderr)
            if process.returncode:
                phase.update(status='checker_failed',exit_code=process.returncode);io.dump(folder/f'{name}.phase.json',phase)
                raise RuntimeError(f'Independent checker failed; see {folder/name}.checker.log')
            receipt=json.loads(rp.read_text());radius=F(receipt['queries'][0]['radius'])
            phase.update(status='target_met' if radius<=tolerance else 'target_not_met',radius=str(radius),receipt=rp.name,candidate=cp.name,verification_wall_seconds=perf_counter()-verify_start,phase_seconds=perf_counter()-began)
            phases.append(phase);io.dump(folder/f'{name}.phase.json',phase)
            if best is None or radius<best[0]:best=(radius,x.copy(),rp,cp)
            if verbose:print(json.dumps(dict(query=i,bond=cap,exact_radius=float(radius),status=phase['status'])),flush=True)
            if radius<=tolerance:break
        if best is None:raise RuntimeError('no construction under supplied budget')
        radius,x,rp,cp=best;warm=x.copy();row=dict(frequency=frequency,status='target_met' if radius<=tolerance else 'target_not_met',radius=str(radius),receipt=str(rp.relative_to(out)),candidate=str(cp.relative_to(out)),phases=phases);results.append(row);io.dump(folder/'RESULT.json',row)
    factor=1 if sys.platform=='darwin' else 1024;parent=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*factor;child=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*factor
    decision=dict(status='target_met' if all(r['status']=='target_met' for r in results) else 'target_not_met',target_radius=str(tolerance),max_radius=str(max(F(r['radius']) for r in results)),queries=results,model_hash=exact.model_hash(spec),sites=n,formal_sector_dimension=comb(n,n//2)**2,global_sector_enumerated=False,determinants_enumerated=0,elapsed_seconds=perf_counter()-start,parent_peak_rss_bytes=parent,child_peak_rss_bytes=child,peak_rss_conservative_sum_bytes=parent+child,learning_scope='per-instance tensor optimization, no pretrained transfer or novelty claim')
    io.dump(out/'DECISION.json',decision);return decision
