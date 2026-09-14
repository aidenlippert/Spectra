"""Stdlib replay of the parity rank profile."""
import json,sys,hashlib
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.discovery.h8_response_rank import valence_states,gram_leakage
from experiments.marginal_determinant_tree import DeterminantOracle
def main():
 d=Path(__file__).resolve().parents[3]; c=json.loads((Path(__file__).parent/'certificate.json').read_text()); h=json.loads((d/'results/marginal_graded_hubbard8/hamiltonian.json').read_text())
 digest=hashlib.sha256(json.dumps(h,sort_keys=True,separators=(',',':')).encode()).hexdigest()
 if digest!='78365286838362d1b591e7afcef0f1e2a085fafc445e961f3b7e19603dbf1c6c' or digest!=c['hamiltonian_canonical_sha256']: raise ValueError('digest')
 p=valence_states(); o=DeterminantOracle(h); g=gram_leakage(o,p); ix=c['parity_even_indices']
 expected=[i for i,s in enumerate(p) if sum(k for k in range(8) if s&(1<<(2*k)))%2==0]
 if ix!=expected or len(ix)!=38 or any(g[i][j] for i in ix for j in ix if i!=j) or min(g[i][i] for i in ix)<2: raise ValueError('profile')
 if any(o.action(s).get(t,F(0)) for s in p for t in p):raise ValueError('Nonzero retained block')
 if c['gamma']!='-381/100' or c['tau']!='-21/5': raise ValueError('threshold')
 delta=F(c['gamma'])-F(c['tau']);minimum=min(g[i][i]/delta+F(c['tau']) for i in ix)
 if delta<=0 or minimum<=0:raise ValueError('Negative restriction failed')
 (Path(__file__).parent/'independent_replay.json').write_text(json.dumps({'status':'verified','stdlib_only':True,'rank_lower_bound':38,'actions':len(o.cache),'minimum_negative_schur_diagonal':str(minimum)},indent=2)+'\n')
if __name__=='__main__': main()
