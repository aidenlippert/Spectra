"""Costed adaptive truncation proposals with independently checked residuals.

The magnitude rule is conventional. The work-weighted rule is a supplied
headroom candidate, not acquired knowledge. Neither changes the tolerance.
"""
from fractions import Fraction as F
from .v7_certificate import Piece,norm_witness, clean


def pruned_taylor(gen,initial,duration,tolerance,*,mode='magnitude',allocation=F(1,4),max_order=24):
 if mode not in ('magnitude','work'):raise ValueError('mode')
 if not 0<allocation<1 or duration<=0 or tolerance<=0 or not 1<=max_order<=24:raise ValueError('budget')
 cs=[clean(initial,gen.n,gen.max_terms)];residual_groups={}; spent=F(0); probes=[];cost=dict(sort_terms=0,group_comparisons=0,rank_column_queries=0)
 slice_budget=tolerance*allocation/F(max_order)
 for m in range(max_order+1):
  raw=gen.apply(cs[-1]);tailfactor=duration**(m+1)/F(m+1)
  for grouping in ('l1','firstfit','weighted'):
   gs,cc=norm_witness(raw,grouping);cost['group_comparisons']+=cc['group_comparisons'];cost['sort_terms']+=cc['sorted_terms']
   bound=spent+tailfactor*sum((F(g['upper']) for g in gs),F(0))
   probes.append(dict(order=m,bound=str(bound),grouping=grouping))
   if bound<=tolerance:
    ws={'jump:0':[],**residual_groups}
    ws[f'residual:0:{m}']=gs
    # Match checker degree trimming when final residual is identically zero.
    if not raw:
     last=max((j for j in range(m) if residual_groups.get(f'residual:0:{j}')),default=0)
     ws={k:v for k,v in ws.items() if k=='jump:0' or int(k.rsplit(':',1)[1])<=last}
    witness=dict(schema='v7-residual-1',integration_basis='power',witnesses=ws,claimed_bound=str(bound),construction_cost=cost)
    return Piece(duration,tuple(cs)),witness,dict(probes=probes,pruning_spent=str(spent),cost=cost,generator=dict(gen.cost),polynomial_entries=sum(map(len,cs)))
  if m==max_order:break
  candidate={p:c/F(m+1) for p,c in raw.items()};timeweight=duration**(m+1)
  costs={p:1 for p in candidate}
  if mode=='work':
   for p in candidate:
    costs[p]=1+len(gen.apply({p:F(1)}));cost['rank_column_queries']+=1
  ranked=sorted(candidate,key=lambda p:(abs(candidate[p])/costs[p],p));cost['sort_terms']+=len(ranked)
  removed={};used=F(0)
  for p in ranked:
   loss=abs(candidate[p])*timeweight
   if used+loss<=slice_budget:
    used+=loss;removed[p]=candidate[p]
  for p in removed:del candidate[p]
  spent+=used
  # R_m=-(m+1) drop_(m+1); sign does not affect norm partition.
  omitted={p:-(m+1)*c for p,c in removed.items()}
  gs,cc=norm_witness(omitted,'l1');cost['sort_terms']+=cc['sorted_terms']
  residual_groups[f'residual:0:{m}']=gs
  cs.append(candidate)
 raise ValueError('pruned Taylor budget without certificate')
