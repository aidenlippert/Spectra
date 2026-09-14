"""Small-Hilbert-space row-sum norm diagnostic; exponential costs explicit."""
from fractions import Fraction as F
from time import perf_counter
import json
from .certificates import _sqrt_interval
from .v7_headroom import ROOT,model,TOL
from .v7_certificate import Generator,norm_witness


def row_bound(op,n):
 if not 1<=n<=6 or len(op)>512:raise ValueError('dense row witness budget')
 dim=2**n;entries={};updates=0
 for p,c in op.items():
  if len(p)!=n:raise ValueError('Pauli width')
  for col in range(dim):
   row=col;phase=0
   for site,char in enumerate(p):
    bit=(col>>(n-1-site))&1
    if char in 'XY':row^=1<<(n-1-site)
    if char=='Y':phase+=1+2*bit
    elif char=='Z':phase+=2*bit
    elif char not in 'IX':raise ValueError('Pauli label')
   re,im=entries.get((row,col),(F(0),F(0)));phase%=4
   if phase==0:re+=c
   elif phase==1:im+=c
   elif phase==2:re-=c
   else:im-=c
   entries[(row,col)]=(re,im);updates+=1
 sums=[F(0)]*dim
 for (row,col),(re,im) in entries.items():
  hi=abs(re)+abs(im) if not re or not im else _sqrt_interval(re*re+im*im,16)[1]
  sums[row]+=hi
 return max(sums),dict(hilbert_dimension=dim,pauli_entry_updates=updates,matrix_entries=len(entries))

def run():
 prior=json.loads((ROOT/'results/v7/strong_baseline.json').read_text())['rows'];rows=[]
 for case in prior:
  n=case['n'];h,o=model(n,case['family']);gamma=F(case['gamma']);T=F(case['T']);g=Generator(h,gamma,n,max_terms=512);c=o
  orders=[case['order']-1] if case['status']=='certified' else list(range(3,8))
  for m in range(max(orders)+1):
   try:r=g.apply(c)
   except ValueError:break
   if m in orders:
    factor=T**(m+1)/F(m+1);start=perf_counter();b,cost=row_bound(r,n);elapsed=perf_counter()-start
    ordinary=[]
    for kind in ('l1','firstfit','weighted'):
     gs,_=norm_witness(r,kind);ordinary.append(sum(F(x['upper']) for x in gs)*factor)
    rows.append(dict(n=n,family=case['family'],gamma=str(gamma),T=str(T),order=m,terms=len(r),row_bound=str(b*factor),group_bound=str(min(ordinary)),
                     candidate_pass=b*factor<=TOL,strict_gain=min(ordinary)>TOL and b*factor<=TOL,cost=cost,elapsed=elapsed))
   c={p:v/F(m+1) for p,v in r.items()}
 result=dict(rows=rows,scope='development norm diagnostic only; 2^n row construction counted; no acquired method or complete-cost gate')
 (ROOT/'results/v7/row_norm_probe.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(rows=len(rows),strict_gains=[{k:r[k] for k in ('n','family','gamma','T','order','terms','elapsed')} for r in rows if r['strict_gain']]),indent=2))
if __name__=='__main__':run()
