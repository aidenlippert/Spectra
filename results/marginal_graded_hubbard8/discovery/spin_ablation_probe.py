"""Propose a physical negative vector without T, then check its exact CAR energy."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
import numpy as np
from scipy.linalg import eigh
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from signed_charge_numeric import prepare,BASE
from experiments.marginal_local_hubbard_block import _actions
from experiments.marginal_hopping_telescope import actions,projected_matrix
from experiments.marginal_spin_telescope import actions as spin_actions
from experiments.marginal_range_two_density import diagonal_value
from experiments.marginal_signed_charge_telescope import local_value as signed_value
from experiments.marginal_quadratic_charge_telescope import local_value as quadratic_value
from experiments.marginal_charge_square_pairs import local_value as square_value
from experiments.marginal_charge_indicator_telescope import local_value as indicator_value


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve()
    cp=folder/'profile_joint_r1_2_certificate.json';rp=folder/'range_two_replay.json';r=json.loads(rp.read_text());c=json.loads(cp.read_text())
    if not r['accepted'] or c['kind']!='hubbard_projector_extension_v12' or r['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted hopping energy certificate required')
    x,_,mats,half,hn,charged,cn,ratio=prepare(c,[])
    hopping_action=actions();gamma=F(c["hopping_telescope"])
    sparse={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};signed={k:F(v) for k,v in c['signed_charge_telescope'].items()};alpha=F(c['penalty']);beta=F(c['joint']['penalty']);ell=F(c['penalized_lower']);best=None
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        d=np.array([float(sparse.get(next(iter(col))&1023,0)-sparse.get(next(iter(col))>>2,0)+signed_value(next(iter(col)),signed)) for col in cols]);scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);h=np.array(projected_matrix(cols,hopping_action),float)/scale[:,None]/scale[None,:];ev,vectors=eigh(a+float(gamma)*h+np.diag(d)+float(alpha)*ph+float(beta)*q,subset_by_index=[0,0])
        if best is None or ev[0]<best[0]:best=(float(ev[0]),key,cols,vectors[:,0])
    _,key,cols,w=best;scale=np.sqrt([sum(a*a for a in col.values()) for col in cols]);w=w/scale;coeff=[round(float(a/max(abs(w)))*10**8) for a in w];vector={s:a*b for a,col in zip(coeff,cols) if a for s,b in col.items()};norm=sum(a*a for a in vector.values())
    # This acceptance calculation contracts fresh CAR actions directly in
    # determinant coordinates; it does not reuse the numerical block matrix.
    local=c['local_window'];u=list(map(F,local['onsite_profile']));hop=list(map(F,local['hopping_profile']));density=list(map(F,local['density_profile']));range2=list(map(F,local['range_two_density_profile']))
    physical=_actions(6,F(local['U']),F(local['t']),u,hop,F(local['V']),density)
    base=F(sum(F(a)*b*vector.get(t,0) for s,a in vector.items() for t,b in physical[s].items()),norm)
    diagonal=F(0)
    for s,a in vector.items():
        value=sparse.get(s&1023,0)-sparse.get(s>>2,0)+diagonal_value(s,range2)
        for fn,name in [(quadratic_value,'quadratic_charge_telescope'),(square_value,'charge_square_pair_telescope'),(indicator_value,'higher_charge_indicator_telescope'),(signed_value,'signed_charge_telescope')]:value+=fn(s,{k:F(v) for k,v in c[name].items()})
        diagonal+=F(a*a,norm)*value
    def fidelity(v,n):return F(sum(a*v.get(s,0) for s,a in vector.items())**2,norm*n)
    ph=fidelity(half,hn);q=ph+ratio*sum((fidelity(v,cn) for v in charged),F(0));without=base+diagonal+alpha*ph+beta*q-ell+gamma*F(sum(a*b*vector.get(t,0) for s,a in vector.items() for t,b in hopping_action[s].items()),norm)
    action=spin_actions({key:F(value) for key,value in c['spin_telescope'].items()});change=F(sum(a*b*vector.get(t,0) for s,a in vector.items() for t,b in action[s].items()),norm)
    if not without<0<=without+change:raise ValueError('Exact ablation does not separate local positivity')
    files={Path(__file__).resolve(),Path(__file__).with_name('signed_charge_numeric.py'),Path(__file__).with_name('joint_profile_numeric.py'),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'sector':list(key),'integer_vector':{str(s):a for s,a in vector.items()},'vector_norm':norm,'without_spin_residual':str(without),'without_spin_residual_float':float(without),'with_spin_residual':str(without+change),'with_spin_residual_float':float(without+change),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Exact negative Rayleigh witness for the submitted local lower after deleting only the spin telescopes. Same vector has nonnegative residual with it. This proves necessity at these fixed certificate coefficients, not strict separation from the fully reoptimized older family.'}
    (folder/'spin_ablation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['integer_vector','source_sha256']}))


if __name__=='__main__':main()
