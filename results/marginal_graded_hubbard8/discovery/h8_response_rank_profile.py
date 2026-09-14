"""Parity profile giving a 38-dimensional exact response-rank obstruction."""
from pathlib import Path
import json, hashlib, sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.discovery.h8_response_rank import valence_states,gram_leakage
from experiments.marginal_determinant_tree import DeterminantOracle
ROOT=Path(__file__).resolve().parents[1]
def main(out):
 h=json.loads((ROOT/'hamiltonian.json').read_text()); o=DeterminantOracle(h); p=valence_states(); g=gram_leakage(o,p)
 parity=[sum(i for i in range(8) if (s>>(2*i))&1)%2 for s in p]; even=[i for i,x in enumerate(parity) if x==0]
 if len(even)!=38 or any(g[i][j] for i in even for j in even if i!=j) or min(g[i][i] for i in even)<2: raise ValueError('parity profile failed')
 c={'schema':'h8-response-rank-profile-v1','hamiltonian_canonical_sha256':hashlib.sha256(json.dumps(h,sort_keys=True,separators=(',',':')).encode()).hexdigest(),'gamma':'-381/100','tau':'-21/5','retained_states':p,'parity_even_indices':even,'parity_sizes':[38,32],'even_diagonal_counts':{str(x):sum(g[i][i]==x for i in even) for x in sorted(set(g[i][i] for i in even))},'rank_lower_bound':38,'exact_offdiag_zero':True,'minimum_even_diagonal':str(min(g[i][i] for i in even)),'scope':'Exact coordinate principal restriction; response corrections of rank at most 37 cannot certify this target.'}
 out=Path(out);out.mkdir(parents=True,exist_ok=True);(out/'certificate.json').write_text(json.dumps(c,indent=2)+'\n')
if __name__=='__main__': main(str(ROOT/'response_rank_profile'))
