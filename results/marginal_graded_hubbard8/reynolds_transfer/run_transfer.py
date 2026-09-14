"""Replay the fixed Reynolds certificate after a symmetry-breaking onsite perturbation."""
from fractions import Fraction as F
import copy, hashlib, json, subprocess, sys
from pathlib import Path
from experiments.marginal_polynomial_metric import replay

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'fixed_reynolds/certificate.json'

def perturb(certificate, epsilon):
    out=copy.deepcopy(certificate)
    terms=out['hamiltonian']
    for mode in (0,1):
        terms.append({'word':[[1,mode],[0,mode]],'coefficient':str(epsilon)})
    return out

def hopping_perturbation(certificate, epsilon):
    out=copy.deepcopy(certificate)
    # Left edge is spatial sites 0--1: modes 0<->2 and 1<->3.
    for a,b in ((0,2),(1,3)):
        out['hamiltonian'].extend([
            {'word':[[1,a],[0,b]],'coefficient':str(epsilon)},
            {'word':[[1,b],[0,a]],'coefficient':str(epsilon)}])
    return out

def main():
    raw=SOURCE.read_bytes(); source=json.loads(raw)
    records=[]
    cases=[('onsite_plus_1e-4','onsite',F(1,10000),perturb(source,F(1,10000))),
           ('onsite_minus_1e-4','onsite',F(-1,10000),perturb(source,F(-1,10000))),
           ('onsite_plus_1e-3','onsite',F(1,1000),perturb(source,F(1,1000))),
           ('onsite_plus_1e-2','onsite',F(1,100),perturb(source,F(1,100))),
           ('hopping_plus_1e-4','hopping',F(1,10000),hopping_perturbation(source,F(1,10000))),
           ('hopping_minus_1e-4','hopping',F(-1,10000),hopping_perturbation(source,F(-1,10000)))]
    records=[{'case':'original','kind':'baseline','epsilon':'0','status':'accepted','error':None,'receipt':replay(source)}]
    for name,kind,epsilon,candidate in cases:
        try:
            receipt=replay(candidate); status='accepted'; error=None
        except ValueError as exc:
            receipt=None; status='rejected'; error=str(exc)
        cert_path=HERE/(name+'.json')
        cert_path.write_text(json.dumps(candidate,indent=2)+'\n')
        record={'case':name,'kind':kind,'epsilon':str(epsilon),'status':status,'error':error,'receipt':receipt,'certificate':str(cert_path)}
        if status=='accepted':
            replay_path=HERE/(name+'_independent_replay.json')
            code="import sys,json;sys.path.insert(0,'.');from pathlib import Path;from experiments.marginal_polynomial_metric import replay;p=Path(sys.argv[1]);Path(sys.argv[2]).write_text(json.dumps(replay(json.loads(p.read_text())),indent=2)+'\\n')"
            subprocess.run([sys.executable,'-S','-c',code,str(cert_path),str(replay_path)],cwd=HERE.parents[2],check=True)
            record['independent_replay']=str(replay_path)
        records.append(record)
    result={'source':str(SOURCE),'source_sha256':hashlib.sha256(raw).hexdigest(),
            'perturbation':'Onsite epsilon*(n_0+n_1) and equal Hermitian left-edge hopping epsilon between (0,2),(1,3); exact rational epsilon',
            'records':records,
            'scope':'Transfer of the unchanged fixed proof under symmetry-breaking onsite and nonmonotone left-edge hopping perturbations; larger rejected amplitudes measure fixed-proof margin only, not physical impossibility; no new energy interval or rediscovery claim.'}
    HERE.mkdir(exist_ok=True); (HERE/'transfer.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'source_sha256':result['source_sha256'],'records':[(x['epsilon'],x['status'],x['error']) for x in records]},indent=2))
if __name__=='__main__': main()
