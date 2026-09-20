"""Certify the best of three declared Hubbard rung couplings.

This is a finite parameter-selection demonstration, not synthesis planning.
The frozen response seed and its construction cost are reported separately.
"""
from copy import deepcopy
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter
from hashlib import sha256
import argparse,json,subprocess,sys,resource
from . import io,exact
from .variational import sweep_solve,enrich
from .solver import screen

def main():
    p=argparse.ArgumentParser();p.add_argument('--seed',type=Path,required=True);p.add_argument('--request',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    raw=a.seed.read_bytes();candidate=json.loads(raw);request=json.loads(a.request.read_text());spec=request['model'];z=2+1.5j
    if candidate['model_hash']!=exact.model_hash(spec):raise ValueError('seed model mismatch')
    seed=io.import_program(candidate['queries'][0]['program']);a.out.mkdir(parents=True);start=perf_counter();rows=[]
    for rung in ['3/4','1','5/4']:
        model=deepcopy(spec)
        for edge in model['edges']:
            if edge[0]%2==0 and edge[1]==edge[0]+1:edge[2]=rung
        x=seed.copy();phases=[]
        for cap in [96,160]:
            if x.bond<cap:x,growth=enrich(model,z,x,cap)
            x,phase=sweep_solve(model,z,x,sweeps=3,target=.0004);phases.append(phase)
            if screen(model,x,z)['radius']<=.0008:break
        folder=a.out/rung.replace('/','_');folder.mkdir();io.dump(folder/'request.json',io.request(model,[z]));io.dump(folder/'candidate.json',io.candidate(model,[z],[x]))
        proc=subprocess.run([sys.executable,'-B','-S','-m','spectra_tensor.exact',str((folder/'candidate.json').resolve()),'--request',str((folder/'request.json').resolve()),'--out',str((folder/'receipt.json').resolve())],capture_output=True,text=True,cwd=Path(__file__).resolve().parent.parent)
        (folder/'checker.log').write_text(proc.stdout+proc.stderr)
        if proc.returncode:raise RuntimeError(proc.stderr)
        receipt=json.loads((folder/'receipt.json').read_text());row=dict(rung=rung,**receipt['queries'][0],phases=phases,verification_seconds=receipt['seconds']);rows.append(row)
        print('inverse',rung,float(F(row['absorption_lower'])),float(F(row['absorption_upper'])),flush=True)
    winner=max(rows,key=lambda r:F(r['absorption_lower']));separated=all(F(winner['absorption_lower'])>F(r['absorption_upper']) for r in rows if r is not winner)
    result=dict(status='certified_best_among_three' if separated else 'intervals_overlap',winner=winner['rung'] if separated else None,objective='maximize -Im b^T(z-H)^-1 b at z=2+1.5i',results=rows,all_response_targets_met=all(F(r['radius'])<=F(1,1000) for r in rows),seed_sha256=sha256(raw).hexdigest(),seed_construction_cost_included=False,seed_cost_location='cold20/query_00: include both failed bond-48 and accepted bond-96 phases',adaptation_and_check_seconds=perf_counter()-start,parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),global_sector_enumerated=False,scope='Three preset finite Hubbard models; not a global design optimum or manufacturing recipe.')
    io.dump(a.out/'DECISION.json',result)
if __name__=='__main__':main()
