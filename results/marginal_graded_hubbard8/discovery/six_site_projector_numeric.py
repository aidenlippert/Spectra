"""Bounded numerical penalty proposal; exact Gram/PSD replay accepts it."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_projector_extendibility import overlap_grams
from experiments.marginal_local_hubbard_block import sector_matrices
BASE=ROOT/'results/marginal_graded_hubbard8'

def main():
    start=time.monotonic()
    vector=json.loads((BASE/'density_transfer/window_ceiling_certificate.json').read_text())['physical_state']
    local=json.loads((BASE/'density_transfer/certificate.json').read_text())['targets'][0]['lower_certificate']
    ceilings=[]
    for m in (2,3,4):
        norm,grams,_=overlap_grams(vector,m,6)
        value=max(float(eigh(np.array(G,dtype=float)/norm,subset_by_index=[len(G)-1,len(G)-1],eigvals_only=True)[0]) for G in grams.values())
        ceiling=F(math.ceil(value*10**9)+1,10**9)
        ceilings.append({'windows':m,'numerical_sum_maximum':value,'rational_ceiling':str(ceiling),'theta':float(ceiling/m)})
    amps={int(s):a for s,a in vector.items()};norm=sum(a*a for a in amps.values())
    sectors=sector_matrices(6,F(local['U']),F(local['t']),[F(x) for x in local['onsite_profile']],
        [F(x) for x in local['hopping_profile']],F(local['V']),[F(x) for x in local['density_profile']])
    matrices=[]
    for key,(_,K,cols) in sectors.items():
        scales=np.sqrt([sum(a*a for a in v.values()) for v in cols])
        A=np.array(K,dtype=float)/scales[:,None]/scales[None,:]
        weights=np.array([sum(a*amps.get(s,0) for s,a in v.items()) for v in cols],dtype=float)/scales/np.sqrt(float(norm))
        matrices.append((key,A,np.outer(weights,weights)))
    def minimum(k):
        return min((float(eigh(A+k*P,subset_by_index=[0,0],eigvals_only=True)[0]),key) for key,A,P in matrices)
    results=[]
    for item in ceilings:
        theta=item['theta']
        opt=minimize_scalar(lambda k:-(minimum(k)[0]-k*theta)/5,bounds=(0,2),method='bounded',options={'maxiter':60,'xatol':1e-8})
        k=F(format(float(opt.x),'.6f'));value,key=minimum(float(k))
        lower=F(math.floor(value*10**7)-1,10**7)
        results.append({**item,'penalty':str(k),'penalized_lower':str(lower),
            'numerical_minimum':value,'active_sector':list(key),
            'proposed_periodic_density':str((lower-k*F(item['rational_ceiling'])/item['windows'])/5)})
    best=max(results,key=lambda r:F(r['proposed_periodic_density']))
    c={'kind':'hubbard_projector_extension_v2','chain_sites':1000000,'target':{'U':'4','t':'1','V':'1/2'},
        'local_window':local,'vector':vector,'windows':best['windows'],'projector_sum_ceiling':best['rational_ceiling'],
        'penalty':best['penalty'],'penalized_lower':best['penalized_lower']}
    out=BASE/'six_site_projector';out.mkdir(exist_ok=True)
    (out/'certificate.json').write_text(json.dumps(c,indent=2)+'\n')
    (out/'numeric_proposal.json').write_text(json.dumps({'proposals':results,'seconds':time.monotonic()-start,
        'scope':'Numerical search only. Reflection Gram normalization retained. Exact all-state overlap and all-Fock positivity replay required.'},indent=2)+'\n')
    print(json.dumps(results),flush=True)

if __name__=='__main__':main()
