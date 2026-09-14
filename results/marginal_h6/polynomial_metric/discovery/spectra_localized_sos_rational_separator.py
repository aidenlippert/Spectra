import json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from experiments.marginal_polynomial_metric import JointPolynomial,check_nonsingleton_separator
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');out=root/'localized_quadratic_sos_obstruction';out.mkdir(exist_ok=True)
c=json.loads((root/'candidate.json').read_text());o=JointPolynomial(c);compiled=json.loads((root/'compiled.json').read_text());proposal=json.loads((root/'localized_quadratic_sos/degree2_proposal.json').read_text());states=proposal['states'];rounded=[round(x*10**9) for x in proposal['dual_weights']];ionic=[any(((s>>(2*i))&3)==3 for i in range(o.sites)) for s in states];count=sum(ionic);rounded[ionic.index(True)]+=10**9-sum(rounded)
def exact_psd(matrix):
 a=[list(row) for row in matrix];previous=1;rank=0;zeros=0;start=time.monotonic()
 for k in range(len(a)):
  pivot=a[k][k]
  if pivot<0:return {'psd':False,'negative_pivot':k,'rank':rank,'seconds':time.monotonic()-start}
  if not pivot:
   if any(a[k][j] for j in range(k+1,len(a))):return {'psd':False,'nonzero_null_row':k,'rank':rank,'seconds':time.monotonic()-start}
   zeros+=1;continue
  for i in range(k+1,len(a)):
   for j in range(i,len(a)):
    num=pivot*a[i][j]-a[i][k]*a[k][j]
    if num%previous:raise ValueError('Nonexact Bareiss quotient')
    a[i][j]=a[j][i]=num//previous
  previous=pivot;rank+=1
 return {'psd':True,'rank':rank,'nullity':zeros,'dimension':len(a),'seconds':time.monotonic()-start}
for mix in [10000,1000,100]:
 denominator=10**9*count*mix;weights=[(mix-1)*count*w+(10**9 if ionic[i] else 0) for i,w in enumerate(rounded)];mom={}
 for s,w in zip(states,weights):
  mask=s
  while True:
   mom[mask]=mom.get(mask,0)+w
   if not mask:break
   mask=(mask-1)&s
 basis=[sum(1<<i for i in inds) for k in range(3) for inds in combinations(range(o.modes),k)];matrix=[[mom.get(i|j,0) for j in basis] for i in basis];psd=exact_psd(matrix);print('mix',mix,psd,flush=True)
 if not psd['psd']:continue
 charge_mom={}
 for state,weight in zip(states,weights):
  charge=sum(((state>>(2*i))&3)==3 for i in range(o.sites))-1;mask=state
  while True:
   charge_mom[mask]=charge_mom.get(mask,0)+weight*charge
   if not mask:break
   mask=(mask-1)&state
 charge_psd=exact_psd([[charge_mom.get(i|j,0) for j in basis] for i in basis]);print('charge',charge_psd,flush=True)
 if not charge_psd['psd']:continue
 functional={'states':states,'weights':[str(F(w,denominator)) for w in weights]}
 try:receipt=check_nonsingleton_separator(o,{m:int(v) for m,v in compiled['numerator']},compiled['receipt']['numerator_scale'],functional)
 except ValueError as exc:print(str(exc),flush=True);continue
 receipt.update(moment_matrix=psd,charge_moment_matrix=charge_psd,mixture_denominator=mix);c['functional']=functional;c['kind']='joint_polynomial_quadratic_sos_separator_v1';(out/'candidate.json').write_text(json.dumps(c,indent=2)+'\n');(out/'discovery_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True);break
