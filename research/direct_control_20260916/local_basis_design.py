"""Small one-particle localization proposal; no many-electron solve."""
import json,time,argparse
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from pyscf import gto,scf,ao2mo
from experiments.marginal_symbolic import canonical,decode


def design(source,out):
 start=time.monotonic();data=json.loads(Path(source).read_text());m=data['modes'];n=m//2
 mol=gto.M(atom=data['geometry'],basis=data['basis'],unit=data['unit'],spin=0,verbose=0)
 mf=scf.RHF(mol).run(conv_tol=1e-12,verbose=0)
 if not mf.converged:raise ValueError('RHF proposal did not converge')
 c=mf.mo_coeff;S=mol.intor_symmetric('int1e_ovlp');vals,vecs=np.linalg.eigh(S)
 target=c.T@((vecs*np.sqrt(vals))@vecs.T)
 h=c.T@mf.get_hcore()@c;eri=ao2mo.restore(1,ao2mo.kernel(mol,c),n).reshape((n,)*4)
 raw={}
 for p in range(m):
  for q in range(m):
   if p%2==q%2:raw[((1,p),(0,q))]=F(float(h[p//2,q//2]))
   for r in range(m):
    for s in range(m):
     v=(eri[p//2,r//2,q//2,s//2]*(p%2==r%2)*(q%2==s%2)-eri[p//2,s//2,q//2,r//2]*(p%2==s%2)*(q%2==r%2))/4
     if v:raw[((1,p),(1,q),(0,s),(0,r))]=F(float(v))
 raw=canonical(raw);given=decode(data['hamiltonian'],m,4)
 # Orbital phases are only a proposal alignment. Original coefficients remain
 # the mathematical input and will be transformed directly, not regenerated.
 best=None
 for mask in range(1<<(n-1)):
  signs=[1]+[-1 if mask>>(i-1)&1 else 1 for i in range(1,n)]
  e=sum(abs(float(raw.get(w,0))*np.prod([signs[p//2] for _,p in w])-float(given.get(w,0))) for w in set(raw)|set(given))
  if best is None or e<best[0]:best=(e,signs)
 target=np.diag(best[1])@target
 result={'target':target.tolist(),'phase_alignment':best[1],'proposal_coefficient_l1_difference':best[0],
         'seconds':time.monotonic()-start,'many_body_discovery':False,'scope':'orbital proposal only; original rational model is authoritative'}
 Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='target'}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('out');a=p.parse_args();design(a.source,a.out)
