"""Exact all-Fock identities explaining the single new quadratic dual row."""
from pathlib import Path
from fractions import Fraction as F
import json
import hashlib
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_joint_family_limit import _profiles
from experiments.marginal_local_hubbard_block import _actions
from experiments.marginal_quadratic_charge_telescope import LABELS,local_value
from experiments.marginal_range_two_density import diagonal_value


def main():
    baseline=_actions(6,F(10,3),1,*_profiles([F(0)]*6)[:2],F(1,2),_profiles([F(0)]*6)[2])
    gradients={}
    for index,name in [(0,'onsite_a'),(1,'onsite_b'),(4,'density_d'),(5,'density_e')]:
        x=[F(0)]*6;x[index]=F(1)
        u,t,v=_profiles(x)
        varied=_actions(6,F(10,3),1,u,t,F(1,2),v)
        values=[]
        for s in range(4096):
            diff={k:varied[s].get(k,F(0))-baseline[s].get(k,F(0)) for k in set(varied[s])|set(baseline[s])}
            if any(v for k,v in diff.items() if k!=s):raise ValueError('Diagonal profile derivative has off-diagonal terms')
            values.append(diff.get(s,F(0)))
        gradients[name]=values
    gradients['free_range_two']=[diagonal_value(s,[1,-1,-1,1]) for s in range(4096)]
    gradients['new_quadratic']=[local_value(s,{'0,3':F(1)}) for s in range(4096)]
    identities={'0,0':{'onsite_a':2,'onsite_b':-2},'0,1':{'density_d':1,'density_e':-1},
      '0,2':{'free_range_two':1},'0,3':{'new_quadratic':1},'1,1':{'onsite_b':2},'1,2':{'density_e':1}}
    for label,coefficients in identities.items():
        if any(local_value(s,{label:F(1)})!=sum(a*gradients[name][s] for name,a in coefficients.items()) for s in range(4096)):
            raise ValueError('Quadratic moment identity fails')
    # Exact row elimination proves that all six physical telescope columns
    # are independent; the new one is outside the five existing directions.
    pivots={}
    for s in range(4096):
        row=[local_value(s,{label:F(1)}) for label in LABELS]
        for j in sorted(pivots):
            a=row[j];row=[v-a*b for v,b in zip(row,pivots[j])]
        j=next((j for j,v in enumerate(row) if v),None)
        if j is not None:
            a=row[j];pivots[j]=[v/a for v in row]
    if len(pivots)!=6:raise ValueError('Unexpected quadratic telescope rank')
    files={Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'determinants_checked':4096,'quadratic_rank':6,'existing_quadratic_rank':5,
      'identities':identities,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
      'scope':'Exact identities on all six-site Fock determinants against freshly reconstructed CAR profile derivatives. Five quadratic constraints are already covered by nearest and free range-two profiles; exactly one additional independent row is needed. This is a finite operator identity, not representability.'}
    out=ROOT/'results/marginal_graded_hubbard8/quadratic_charge_telescope/quadratic_profile_identities.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__=='__main__':main()
