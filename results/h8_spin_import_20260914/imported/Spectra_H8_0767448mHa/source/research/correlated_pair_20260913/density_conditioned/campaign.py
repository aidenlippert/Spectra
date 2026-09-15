import json
from pathlib import Path
from fractions import Fraction as F
from maxflow import fixture_terms,feasible,shifted_particle_hole
from closure import factor,closure,check_all_terms

def run():
 out=Path('results/correlated_pair_20260913/density_conditioned');out.mkdir(parents=True,exist_ok=True)
 summary={}
 for name in ('h6','h8'):
  d=json.loads(Path(f'results/molecular_collective_20260913/campaign/{name}/fixture.json').read_text());c,e,U,mu=shifted_particle_hole(d)
  edges=[(i,j,-v) for (i,j),v in U.items() if v<0]
  summary[name]={'particle_hole_mu':str(mu),'positive_capacities':list(map(str,e.values())),'negative_edges':len(edges),'flow':feasible(list(e.values()),edges),'pair_terms':len(U)}
 f=factor(((1,0),),(1,1),{((1,2),):F(1)},F(1),F(1),F(1),{((1,3),(1,4),(0,3)):F(1)})
 summary['resonant_closure']=check_all_terms(closure(f))
 (out/'diagnostic.json').write_text(json.dumps(summary,indent=2,default=str))
if __name__=='__main__':run()
