"""Finite sparse Krylov upper proposals; original rational Rayleigh gate accepts."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.discovery.singlet_energy_replay import apply,gram
from experiments.marginal_determinant_tree import DeterminantOracle
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'singlet_upper_krylov';out.mkdir(exist_ok=True)
    c=json.loads((ROOT/'singlet_valence_embedding/certificate.json').read_text())
    h=json.loads((ROOT/'hamiltonian.json').read_text());o=DeterminantOracle(h)
    V=[{int(s):a for s,a in v.items()} for v in c['basis']]
    K=[];powers=V;results=[]
    for degree in range(5):
        try:
            K+=powers;images=[apply(o,k) for k in K]
            G=np.array(gram(K,K),float);A=np.array(gram(K,images),float)
            vals,T=eigh(G);keep=vals>max(vals)*1e-12;T=T[:,keep]/np.sqrt(vals[keep])
            energies,ev=eigh(T.T@A@T);weights=T@ev[:,0]
            v={}
            for a,k in zip(weights,K):
                for s,x in k.items():v[s]=v.get(s,0.)+a*float(x)
            ints={s:round(x*1e12) for s,x in v.items()};states=sorted(s for s,x in ints.items() if x)
            witness={'states':states,'amplitudes':[ints[s] for s in states]}
            upper=o.upper(witness)
            row={'degree':degree,'directions':len(K),'numerical_rank':sum(keep).item(),'upper':str(upper),
                 'upper_float':float(upper),'support':len(states),'unique_action_states':len(o.cache),'referenced_determinants':o.referenced_state_count()}
            results.append(row);(out/f'upper_degree{degree}.json').write_text(json.dumps(witness)+'\n')
            print(row,flush=True);powers=images[-len(V):]
        except ValueError as error:
            results.append({'degree':degree,'refused':str(error),'unique_action_states':len(o.cache)});break
        finally:(out/'receipt.json').write_text(json.dumps({'results':results,'scope':'Finite sparse variational upper proposals, exact original Rayleigh acceptance; no full-sector eigensolve.'},indent=2)+'\n')
if __name__=='__main__':main()
