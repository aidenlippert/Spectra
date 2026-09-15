"""Replayable request gate using the actual v6-selected affine ARX family.

Only the constant-input step record is read. A future Q1=75 intervention is
specified directly; no sine response outcomes are opened or consumed.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import csv, json
import numpy as np
from .v6_headroom import Trace, fit
from .v7_prediction_requests import request_prediction, _feature
from .v7_identifiability import verify_compatible_models, verify_row_certificate

ROOT=Path(__file__).resolve().parents[1]

def run():
 path=ROOT/'data/v6/tclab_step_test.csv'
 rows=list(csv.DictReader(path.open()))
 ys=[F(r['T1']) for r in rows]
 us=[[F(r['Q1']),F(r['Q2'])] for r in rows[:-1]]
 def trace(a,b): return Trace(np.array(us[a:b],dtype=float),np.array(ys[a:b+1],dtype=float))
 model=fit([trace(0,540)],[trace(540,720)],'stateful_arx')
 p=model.p; k=2
 # The exact rational model below defines a candidate. It need not be exactly
 # equal to the floating fit; its training compatibility is checked directly.
 ay=[F(str(x)) for x in model.coef[:p]]
 by=[F(str(model.norm.ys*model.coef[p+i*k+j]/model.norm.us[j])) for i in range(p) for j in range(k)]
 ym=F(str(model.norm.ym)); um=[F(str(x)) for x in model.norm.um]
 d=ym+F(str(model.norm.ys*model.coef[-1]))-ym*sum(ay)-sum(by[i*k+j]*um[j] for i in range(p) for j in range(k))
 theta=ay+by+[d]
 X=[_feature(ys[:t+1],us[:t],us[t],p,k) for t in range(p-1,len(us))]
 y=ys[p:]
 eta=max(abs(sum(a*b for a,b in zip(row,theta))-v) for row,v in zip(X,y))
 training=dict(X=X,y=y,theta0=theta,order=p,inputs=k)
 witness=[F(0)]*len(theta);witness[p]=1;witness[-1]=-50
 request=request_prediction(training,ys,us,[F(75),F(0)],F(1,10),eta,witness=witness,horizon=16)
 w=request['witness']
 valid=verify_compatible_models(X,y,w['query'],w['theta_plus'],w['theta_minus'],eta,F(2,5))
 if not valid or request['status']!='no_guarantee': raise AssertionError('ambiguity gate failed')
 # A particular historical row query is supported despite nonunique theta.
 # This demonstrates informativity for a query, not every unchanged-input future.
 weights=[F(1)]+[F(0)]*(len(X)-1)
 supported=verify_row_certificate(X,X[0],weights)
 if not supported: raise AssertionError('rowspace demonstration failed')
 residual_plus=max(abs(sum(a*b for a,b in zip(row,w['theta_plus']))-v) for row,v in zip(X,y))
 residual_minus=max(abs(sum(a*b for a,b in zip(row,w['theta_minus']))-v) for row,v in zip(X,y))
 output=dict(status='request_refused_with_verified_compatible_models',
  source=str(path.relative_to(ROOT)),source_sha256=sha256(path.read_bytes()).hexdigest(),
  selection='existing v6 stateful_arx fitter, step[0:540] fit and step[540:720] selection',
  admitted_class='unconstrained affine ARX with free intercept and lag-major Q1,Q2 coefficients',
  selected_order=p, training_rows=len(X), feature_columns=len(theta),
  eta=str(eta), eta_float=float(eta), eta_semantics='smallest training residual enclosure about proposed rational model; not a certified hardware noise bound',
  requested_Q1=75,requested_Q2=0,requested_horizon=16,tolerance='1/10',
  prediction_disagreement=str(abs(sum(F(a)*(b-c) for a,b,c in zip(w['query'],w['theta_plus'],w['theta_minus'])))),
  maximum_training_residual_plus=str(residual_plus),maximum_training_residual_minus=str(residual_minus),
  nullspace_witness=list(map(str,witness)), theta_plus=list(map(str,w['theta_plus'])),theta_minus=list(map(str,w['theta_minus'])),query=w['query'],
  exact_witness_verified=valid, future_response_outputs_read=False,
  in_rowspace_example=dict(row=0,center=str(y[0]),radius=str(eta),verified=True),
  limitation='Certification encloses model predictions conditional on this admitted class and uncertainty set; no new physical validation or compounding result.')
 return output

if __name__=='__main__':
 result=run();path=ROOT/'results/v7/prediction_requests.json';path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k] for k in ('status','selected_order','training_rows','eta_float','prediction_disagreement','exact_witness_verified')},indent=2))
