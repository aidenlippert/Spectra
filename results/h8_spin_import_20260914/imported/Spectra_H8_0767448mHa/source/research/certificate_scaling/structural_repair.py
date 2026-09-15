"""Sparse symmetry-preserving ideal repair of newly discovered fixed factors."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,time,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from experiments.marginal_symbolic import verify,decode,encode,expand_squares,add,scale,mono,number_shift,product,multiplier_basis
from research.certificate_scaling.adaptive_block_discovery import partition
from research.certificate_scaling.direct_sparse_discovery import sparse_columns

def repair(c):
 start=time.monotonic();before=verify(c);m=c['modes'];n=c['particles'];h=decode(c['hamiltonian'],m,4)
 squares,_=expand_squares(c['blocks'],c['denominator'],m,c.get('operator_degree',3));target=add(h,scale(squares,-1))
 _,sig,_=partition(h,m,'quadratic',True);zero=sig(())
 basis=[p for p in multiplier_basis(m,max_body=2) if all(sig(w)==zero for w in p)]
 polys=[mono(())]+[product(number_shift(m,n),q) for q in basis]
 rows=sorted(set(target)|{w for p in polys for w in p},key=lambda w:(len(w),w))
 rows=[w for w in rows if tuple(i for c,i in w if c)<=tuple(i for c,i in w if not c)]
 lookup={w:i for i,w in enumerate(rows)};A=sparse_columns([{w:v for w,v in p.items() if w in lookup} for p in polys],lookup)
 weights=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in rows])
 eq=sparse.hstack([A,sparse.eye(len(rows)),-sparse.eye(len(rows))],format='csc')
 r=linprog(np.r_[-1.,np.zeros(len(basis)),weights,weights],A_eq=eq,b_eq=[float(target.get(w,0)) for w in rows],
  bounds=[(None,None)]*A.shape[1]+[(0,None)]*(2*len(rows)),method='highs',options={'time_limit':60.})
 if not r.success:raise RuntimeError(r.message)
 out=dict(c);out['b']=str(F(round(r.x[0]*10**12),10**12));out['number_multiplier']=encode(add(*(scale(q,F(round(v*10**12),10**12)) for q,v in zip(basis,r.x[1:A.shape[1]]))))
 after=verify(out);chosen=out if F(after['lower'])>=F(before['lower']) else c
 return chosen,{'before':before,'after':after,'retained_improvement':chosen is out,'wall_seconds':time.monotonic()-start,'coefficient_rows':len(rows),'free_multiplier_count':len(basis),'LP_iterations':r.nit}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();c,r=repair(json.loads(a.input.read_text()));a.out.mkdir(parents=True,exist_ok=True)
 (a.out/'certificate.json').write_text(json.dumps(c,separators=(',',':'))+'\n');(a.out/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'before':r['before']['lower_float'],'after':r['after']['lower_float'],'seconds':r['wall_seconds']}))
