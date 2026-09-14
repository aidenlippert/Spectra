"""Bounded full-spectrum soft-min annealing; untrusted proposal only."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,sys,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from joint_profile_numeric import build,profiles
from experiments.marginal_charged_projectors import charged_vectors
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value
from signed_charge_numeric import prepare
from experiments.marginal_hopping_telescope import actions,projected_matrix
from experiments.marginal_spin_telescope import LABELS,actions as spin_actions
from experiments.marginal_spectator_hopping import LABELS as SPECTATORS,actions as spectator_actions
from experiments.marginal_pair_transfer import LABELS as PAIRS,actions as pair_actions
BASE=ROOT/'results/marginal_graded_hubbard8'



def main():
    parser=argparse.ArgumentParser();parser.add_argument('seed',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();c=json.loads(args.seed.read_text())
    indices=[1,2,4,5,6,7,31,11,20]
    candidates_path=BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json';candidates=json.loads(candidates_path.read_text())['candidates']
    shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in indices]
    x,physical,mats,*_=prepare(c,shapes)
    action=[actions()]+[spin_actions({label:1}) for label in LABELS]+[spectator_actions({label:1}) for label in SPECTATORS]+[pair_actions({label:1}) for label in PAIRS];hop={}
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);hop[key]=np.array([projected_matrix(cols,item) for item in action],float)/scale[:,None]/scale[None,:]
    theta_h=F(c['projector_sum_ceiling'])/c['windows'];theta_j=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
    diagonal={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};taus=[]
    for shape in shapes:
        ratios={diagonal.get(s,F(0))/v for s,v in shape.items()}
        if len(ratios)!=1:raise ValueError('Seed sparse correction is outside the fixed span')
        taus.append(float(ratios.pop()))
    y0=np.array(list(map(float,x[2:4]))+[float(F(c['penalty'])),float(F(c['joint']['penalty']))]+taus+[float(F(c.get('signed_charge_telescope',{}).get(key,0))) for key in PATTERNS]+[float(F(c.get('hopping_telescope',0)))]+[float(F(c.get('spin_telescope',{}).get(label,0))) for label in LABELS]+[float(F(c.get('spectator_hopping',{}).get(label,0))) for label in SPECTATORS]+[float(F(c.get('pair_transfer',{}).get(label,0))) for label in PAIRS]+[float(F(c['penalized_lower']))])
    assert len(y0)==89
    history_path=args.seed.parent/'thermal_history.json'
    history=json.loads(history_path.read_text());point=np.array(history['last_optimizer_point']);tau=history['stages'][-1]['temperature']
    assert len(point)==88 and tau>0
    spectra=[];minimum=float('inf')
    for key,a,ds,ph,q,ts,*_ in mats:
        matrix=a+np.einsum('i,ijk->jk',point[:2]-np.array(x[2:4],float),ds)+point[2]*ph+point[3]*q+np.diag(point[4:65]@ts)+np.einsum('i,ijk->jk',point[65:88],hop[key])
        ev,vectors=eigh(matrix);minimum=min(minimum,float(ev[0]));spectra.append((key,ev,vectors,ds,ph,q,ts))
    partition=sum(float(np.exp(-(ev-minimum)/tau).sum()) for key,ev,v,ds,ph,q,ts in spectra)
    atoms=[];mean=np.zeros(88);omitted=0.
    for key,ev,vectors,ds,ph,q,ts in spectra:
        weights=np.exp(-(ev-minimum)/tau)/partition
        rho=(vectors*weights)@vectors.T
        mean+=np.r_[[np.sum(rho*d) for d in ds],np.sum(rho*ph),np.sum(rho*q),ts@np.diag(rho),np.einsum('ij,kij->k',rho,hop[key])]
        for i,weight in enumerate(weights):
            if weight>1e-12:atoms.append({'sector':key,'dual_weight':float(weight)/5,'coordinates':vectors[:,i].tolist()})
            else:omitted+=float(weight)
    residual=mean.copy();residual[2]-=float(theta_h);residual[3]-=float(theta_j)
    result={'accepted':False,'atoms':atoms,'temperature':tau,'omitted_weight':omitted,'evaluation_point':point.tolist(),'moment_residuals':residual.tolist(),'largest_residuals':sorted(enumerate(residual),key=lambda pair:abs(pair[1]),reverse=True)[:12],'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),args.seed,history_path,Path(__file__).with_name('signed_charge_numeric.py')]},'scope':'Untrusted Gibbs eigenvectors at the final unrounded optimizer point. Residuals2,3 are fidelity minus ceiling and need only be nonpositive; other84-coordinate residuals must vanish for the unrestricted family. Exact physical mixture acceptance remains separate.'}
    (args.output/'dual_atoms.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'atoms':len(atoms),'temperature':tau,'largest_residuals':result['largest_residuals']}))


if __name__=='__main__':main()
