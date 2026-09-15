"""Separate, versioned extension for exact quadratic Pauli norm witnesses.

No learned proposal is trusted. This module is an experimental checker extension
before acquisition; it does not change the frozen original residual checker.
"""
from collections import Counter
from fractions import Fraction as F
from .v7_certificate import residual_records, norm_witness, proof_fraction, rational
from .v7_quadratic_probe import square
from .certificates import _anti,_sqrt_interval


def check_groups(op,groups,cost):
 if not isinstance(groups,list) or len(groups)>len(op):raise ValueError('group count')
 seen=Counter();total=F(0)
 for g in groups:
  labels=g['labels']
  if not isinstance(labels,list) or not labels or len(labels)>len(op):raise ValueError('group labels')
  for i,p in enumerate(labels):
   if p not in op:raise ValueError('unknown term')
   for q in labels[i+1:]:
    cost['group_pairs']+=1
    if not _anti(p,q):raise ValueError('invalid group')
  hi=proof_fraction(g['upper']);sq=sum((op[p]**2 for p in labels),F(0));cost['squares']+=len(labels)+1
  if hi<0 or hi*hi<sq:raise ValueError('invalid group bound')
  seen.update(labels);total+=hi
 if seen!=Counter(op.keys()):raise ValueError('group coverage')
 return total


def check(gen,initial,pieces,witness,tolerance,*,expected_time):
 cost=dict(group_pairs=0,square_pairs=0,squares=0)
 try:
  tol=rational(tolerance);T=rational(expected_time)
  if tol<0 or T<=0:raise ValueError('invalid request')
  if witness['schema']!='v7-quadratic-1':raise ValueError('schema')
  records=residual_records(gen,initial,pieces,witness['integration_basis'])
  if sum((rational(p.duration) for p in pieces),F(0))!=T:raise ValueError('horizon mismatch')
  ws=witness['witnesses']
  if set(ws)!={key for key,_,_ in records}:raise ValueError('record coverage')
  bound=F(0)
  for key,op,factor in records:
   proof=ws[key]
   if proof['kind']=='groups':b=check_groups(op,proof['groups'],cost)
   elif proof['kind']=='square':
    sq,pairs=square(op,cap=gen.max_terms);cost['square_pairs']+=pairs
    # Products recomputed from actual residual, never accepted as supplied.
    square_bound=check_groups(sq,proof['groups'],cost)
    b=proof_fraction(proof['upper'])
    if b<0 or b*b<square_bound:raise ValueError('invalid square root')
   else:raise ValueError('norm proof kind')
   bound+=factor*b
  if bound!=proof_fraction(witness['claimed_bound']):raise ValueError('bound mismatch')
  return dict(status='certified' if bound<=tol else 'over_tolerance',bound=str(bound),cost=cost,generator=dict(gen.cost),time=str(T))
 except (ValueError,KeyError,TypeError,IndexError,AttributeError,ZeroDivisionError) as exc:
  return dict(status='rejected',reason=str(exc),cost=cost,generator=dict(gen.cost))


def adaptive_with_square(gen,initial,duration,tolerance,*,trigger_ratio=F(3,2),term_cap=40,max_order=24):
 if not 0<=max_order<=24 or not 1<=term_cap<=512 or trigger_ratio<=1:raise ValueError('candidate limits')
 cs=[dict(initial)];probes=[];quad_tests=0
 for m in range(max_order+1):
  raw=gen.apply(cs[-1]);factor=duration**(m+1)/F(m+1);opts=[]
  for grouping in ('l1','firstfit','weighted'):
   gs,c=norm_witness(raw,grouping);b=sum((F(g['upper']) for g in gs),F(0))*factor
   opts.append((b,dict(kind='groups',groups=gs)));probes.append(dict(order=m,kind=grouping,bound=str(b)))
   if b<=tolerance:break
  bound,proof=min(opts,key=lambda z:z[0])
  if tolerance<bound<=trigger_ratio*tolerance and len(raw)<=term_cap:
   squared,pairs=square(raw,cap=gen.max_terms);quad_tests+=1
   for grouping in ('firstfit','weighted'):
    gs,c=norm_witness(squared,grouping);b2=sum((F(g['upper']) for g in gs),F(0));hi=_sqrt_interval(b2,16)[1];b=hi*factor
    probes.append(dict(order=m,kind='square_'+grouping,bound=str(b),pairs=pairs))
    if b<bound:bound,proof=b,dict(kind='square',groups=gs,upper=str(hi))
    if bound<=tolerance:break
  if bound<=tolerance:
   ws={'jump:0':dict(kind='groups',groups=[])}
   ws.update({f'residual:0:{j}':dict(kind='groups',groups=[]) for j in range(m+1)})
   if raw:ws[f'residual:0:{m}']=proof
   else:ws={'jump:0':dict(kind='groups',groups=[]),'residual:0:0':dict(kind='groups',groups=[])}
   from .v7_certificate import Piece
   return Piece(duration,tuple(cs)),dict(schema='v7-quadratic-1',integration_basis='power',witnesses=ws,claimed_bound=str(bound)),dict(probes=probes,quadratic_tests=quad_tests,generator=dict(gen.cost))
  if m<max_order:cs.append({p:c/F(m+1) for p,c in raw.items()})
 raise ValueError('order budget without certificate')
