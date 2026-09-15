"""Exact reuse opportunity diagnostic; no timing/learned headroom claim."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
from .v7_headroom import model
from .v7_certificate import Generator
ROOT=Path(__file__).resolve().parents[1]

def trajectory(h,gamma,n,initial,B,degree):
 g=Generator(h,gamma,n,512);c=dict(initial);coeffs=[];raw=[];ops=[]
 for k in range(degree+1):
  coeffs.append(c);before=g.cost['coefficient_multiply_adds'];full=g.apply(c)
  raw.append(full);ops.append(g.cost['coefficient_multiply_adds']-before)
  c={p:v/F(k+1) for p,v in full.items() if p in B}
 return coeffs,raw,ops

def run():
 rows=[]
 chosen={(3,'xxz','0','1/5'),(3,'xxz','2','1/2'),(4,'mixed','2','1/2'),(6,'xxz','2','1/5')}
 for item in json.loads((ROOT/'results/v9/headroom.json').read_text())['rows']:
  key=(item['n'],item['family'],item['gamma'],item['time'])
  if item['arm']!='galerkin' or key not in chosen:continue
  h,o=model(item['n'],item['family']);B=set(o);before=None;comparisons=[];all_ops=0;savable=0;entries=0;delta_ops=0;merge_visits=0
  for r in item['cost']['rounds']:
   cs,raw,ops=trajectory(h,F(item['gamma']),item['n'],o,B,r['order']);all_ops+=sum(ops);entries=max(entries,sum(map(len,raw)))
   if before is not None:
    oldcs,oldraw,oldB=before;newlabels=B-oldB
    first=next((k for k,op in enumerate(oldraw) if set(op)&newlabels),len(oldraw))
    reusable=min(first+1,len(oldraw),len(raw))
    for k in range(reusable):
     if cs[k]!=oldcs[k] or raw[k]!=oldraw[k]:raise AssertionError('prefix lemma violated')
    saved=sum(ops[:reusable]);savable+=saved
    dg=Generator(h,F(item['gamma']),item['n'],512);diff_work=0;merge_work=0
    for k,c in enumerate(cs):
     if k>=len(oldcs):diff_work+=ops[k];continue
     old=oldcs[k];labels=set(c)|set(old)
     difference={p:c.get(p,F(0))-old.get(p,F(0)) for p in labels};difference={p:v for p,v in difference.items() if v}
     delta_raw=dg.apply(difference)
     reconstructed=dict(oldraw[k])
     for p,v in delta_raw.items():
      reconstructed[p]=reconstructed.get(p,F(0))+v
      if not reconstructed[p]:del reconstructed[p]
     if reconstructed!=raw[k]:raise AssertionError('difference identity failed')
     merge_work+=len(labels)+len(oldraw[k])+len(delta_raw)
    diff_work+=dg.cost['coefficient_multiply_adds'];delta_ops+=diff_work;merge_visits+=merge_work
    comparisons.append(dict(old_degree=len(oldcs)-1,new_degree=len(cs)-1,first_changed_projection=first,reusable_generator_actions=reusable,savable_multiply_adds=saved,difference_multiply_adds=diff_work,additional_merge_visits=merge_work))
   else:delta_ops+=sum(ops)
   before=(cs,raw,set(B));B.update(r.get('added_labels',[]))
  rows.append(dict(n=key[0],family=key[1],gamma=key[2],time=key[3],rounds=len(item['cost']['rounds']),generator_multiply_adds=all_ops,theoretically_reusable_multiply_adds=savable,difference_generator_multiply_adds=delta_ops,additional_merge_visits=merge_visits,peak_cached_raw_entries=entries,comparisons=comparisons))
 out=dict(schema='v9-reuse-opportunity-1',rows=rows,claim='exact coefficient reuse validated; no implementation cost win or acquisition established',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (ROOT/'results/v9/reuse_probe.json').write_text(json.dumps(out,indent=2)+'\n')
 for r in rows:print(r['n'],r['family'],r['gamma'],r['time'],r['theoretically_reusable_multiply_adds'],r['generator_multiply_adds'],'rawentries',r['peak_cached_raw_entries'])
if __name__=='__main__':run()
