"""Stdlib replay for h8_response_rank certificate (no NumPy)."""
import hashlib,json,sys
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_schur_transfer import ldl_pivots
EXPECTED='78365286838362d1b591e7afcef0f1e2a085fafc445e961f3b7e19603dbf1c6c'
def states():
 return [sum(1<<(2*i) for i in u)|sum(1<<(2*i+1) for i in range(8) if i not in u) for u in combinations(range(8),4)]
def replay(path):
 c=json.loads(Path(path).read_text()); h=json.loads((ROOT/'results/marginal_graded_hubbard8/hamiltonian.json').read_text())
 if hashlib.sha256(json.dumps(h,sort_keys=True,separators=(',',':')).encode()).hexdigest()!=EXPECTED: raise ValueError('Hamiltonian digest mismatch')
 if c.get('gamma')!='-413/100': raise ValueError('gamma mismatch')
 p=states()
 if c.get('retained_states')!=p or len(p)!=70: raise ValueError('valence state mismatch')
 o=DeterminantOracle(h); retained=set(p); gram=[]
 for s in p:
  a={t:v for t,v in o.action(s).items() if t not in retained}
  gram.append([sum(x*b.get(t,F(0)) for t,x in a.items()) for b in [{t:v for t,v in o.action(q).items() if t not in retained} for q in p]])
 if any(o.action(s).get(t,F(0)) for i,s in enumerate(p) for t in p): raise ValueError('H_PP is not zero')
 if [[str(x) for x in r] for r in gram]!=c.get('gram'): raise ValueError('Gram mismatch')
 if type(c.get('tests')) is not list or not 1<=len(c['tests'])<=10: raise ValueError('Bounded nonempty target list required')
 rows=[]
 for row in c.get('tests',[]):
  tau=F(row['tau']); den=F(-413,100)-tau
  if den<=0: raise ValueError('invalid tau')
  ix=row['principal_indices']
  if len(ix)<33 or len(set(ix))!=len(ix) or any(type(i) is not int or not 0<=i<70 for i in ix): raise ValueError('invalid indices')
  s=[[-tau*(i==j)-gram[i][j]/den for j in range(70)] for i in range(70)]
  piv=ldl_pivots([[-s[i][j] for j in ix] for i in ix])
  if piv is None: raise ValueError('claimed restriction is not negative definite')
  rows.append({'tau':row['tau'],'exact_ldl_positive':True,'dimension':len(ix)})
 out=Path(path).with_name('independent_replay.json'); out.write_text(json.dumps({'status':'verified','stdlib_only':True,'tests':rows,'actions':len(o.cache)},indent=2)+'\n'); return rows
if __name__=='__main__': print(replay(sys.argv[1] if len(sys.argv)>1 else str(Path(__file__).with_name('certificate.json'))))
