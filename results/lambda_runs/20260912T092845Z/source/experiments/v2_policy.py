"""Exact finite Bayesian experiment certificates."""
from fractions import Fraction
from math import gcd
from functools import reduce
V=Fraction(1,2); MODELS=tuple((a,b,s) for a in (0,1) for b in (0,1) for s in (-1,1)); ACTIONS=tuple((i,j) for i in (0,1) for j in (0,1))
def _norm(w):
 g=reduce(gcd,(abs(x) for x in w if x),0) or 1; return tuple(x//g for x in w)
def _prior(k): return tuple(1 if (k not in ('A','AB') or a==0) and (k not in ('B','AB') or b==0) else 0 for a,b,s in MODELS)
def _num(m,a,y):
 a0,b0,s=m; q=3 if (a0==a[0] and b0==a[1]) else 2; return q if (s==1)==bool(y) else 4-q
def _child(w,a,y): return _norm(tuple(x*_num(m,a,y) for x,m in zip(w,MODELS)))
def _risk(w):
 p=sum(x for x,m in zip(w,MODELS) if m[2]==1); n=sum(x for x,m in zip(w,MODELS) if m[2]==-1); return Fraction(min(p,n),sum(w)) if sum(w) else Fraction(0)
def certificate(known='none',horizon=1):
 if type(horizon) is not int or known not in ('none','A','B','AB') or not 0<=horizon<=5: raise ValueError('unsupported policy task')
 w0=_prior(known); memo={}; nodes={}
 def solve(w,h):
  key=(w,h)
  if key in memo:return memo[key]
  if h==0 or _risk(w)==0:
   v=_risk(w); memo[key]=v; nodes[key]={'weights':list(w),'horizon':h,'value':str(v),'terminal':True,'actions':[]}; return v
  acts=[]
  for a in ACTIONS:
   bs=[]; val=Fraction(0); total=sum(w)
   for y in (0,1):
    z=_child(w,a,y); p=Fraction(sum(x*_num(m,a,y) for x,m in zip(w,MODELS)),4*total); val+=p*solve(z,h-1); bs.append({'outcome':y,'probability':str(p),'weights':list(z)})
   acts.append({'action':list(a),'value':str(val),'branches':bs})
  v=min(Fraction(x['value']) for x in acts); memo[key]=v; nodes[key]={'weights':list(w),'horizon':h,'value':str(v),'terminal':False,'actions':acts}; return v
 v=solve(w0,horizon); return {'known':known,'horizon':horizon,'models':[list(m) for m in MODELS],'actions':[list(a) for a in ACTIONS],'value':str(v),'nodes':list(nodes.values())}
def verify(c):
 try:
  k=c['known']; h=c['horizon']
  if type(h) is not int: return False
  if k not in ('none','A','B','AB') or not 0<=h<=5 or c['models']!=[list(m) for m in MODELS] or c['actions']!=[list(a) for a in ACTIONS]: return False
  raw=c['nodes']
  if not isinstance(raw,list) or len(raw)>=100000:return False
  table={}
  for n in raw:
   if type(n.get('horizon')) is not int or not 0<=n['horizon']<=h:return False
   ww=n.get('weights')
   if not isinstance(ww,list) or len(ww)!=8 or any(type(x) is not int or x<0 or x.bit_length()>16 for x in ww) or sum(ww)<=0 or tuple(ww)!=_norm(tuple(ww)):return False
   key=(tuple(ww),n['horizon'])
   if key in table:return False
   table[key]=n
  visited=set(); checked=set()
  def visit(w,t):
   if (w,t) in checked:return True
   visited.add((w,t))
   n=table.get((w,t))
   if n is None or sum(w)<=0:return False
   if t==0 or _risk(w)==0:return n['terminal'] is True and n['actions']==[] and n['value']==str(_risk(w))
   if n['terminal'] is not False or len(n['actions'])!=4:return False
   vals=[]; seen=set()
   for act in n['actions']:
    a=tuple(act['action'])
    if a not in ACTIONS or a in seen or len(act['branches'])!=2:return False
    seen.add(a); val=Fraction(0)
    for y,b in enumerate(act['branches']):
     if b['outcome']!=y:return False
     z=tuple(b['weights']); p=Fraction(sum(x*_num(m,a,y) for x,m in zip(w,MODELS)),4*sum(w))
     if b['probability']!=str(p) or z!=_child(w,a,y) or not visit(z,t-1):return False
     val+=p*Fraction(table[(z,t-1)]['value'])
    if act['value']!=str(val):return False
    vals.append(val)
   valid=seen==set(ACTIONS) and n['value']==str(min(vals))
   if valid:checked.add((w,t))
   return valid
  w=_prior(k); return visit(w,h) and c['value']==table[(w,h)]['value'] and visited==set(table)
 except (AttributeError,KeyError,TypeError,ValueError,ZeroDivisionError): return False
def optimal_risk(known='none',horizon=1): return Fraction(certificate(known,horizon)['value'])
