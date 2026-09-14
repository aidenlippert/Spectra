"""H-only orbital graph/tail/quotient counts, without PSD solve or Fock basis."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import decode,adj,canonical,hermitian
from research.certificate_scaling.adaptive_block_discovery import partition

def graph(h,modes):
    # Spatial orbitals, not spin orbitals. Every term induces a support clique.
    n=modes//2;g=[set() for _ in range(n)];cost={}
    for w in h:
        for i,j in combinations(sorted({i//2 for c,i in w}),2):
            g[i].add(j);g[j].add(i);cost[i,j]=cost.get((i,j),F(0))+abs(h[w])
    edge=sum(map(len,g))//2;complete=edge==n*(n-1)//2;active=set(range(n));width=0
    while active:
        def score(i):
            nbr=g[i]&active;return (sum(v not in g[u] for u,v in combinations(nbr,2)),len(nbr),i)
        i=min(active,key=score);nbr=g[i]&active;width=max(width,len(nbr));active.remove(i)
        for u,v in combinations(nbr,2):g[u].add(v);g[v].add(u)
    return {'vertices':n,'edges':edge,'complete':complete,'minimum_edge_deletion_l1':str(min(cost.values())) if cost else '0','treewidth_exact':n-1 if complete else None,'treewidth_minfill_upper':width}

def trim(h,budget):
    seen=set();groups=[]
    for w in h:
        if w in seen:continue
        related=set(canonical(adj({w:F(1)})));related.add(w);seen.update(related)
        groups.append((sum(abs(h.get(u,F(0))) for u in related),related))
    out=dict(h);used=F(0)
    for norm,ws in sorted(groups,key=lambda x:x[0]):
        if used+norm>budget:break
        for w in ws:out.pop(w,None)
        used+=norm
    assert hermitian(out);return out,used

def main():
    rows=[]
    for root in ['active_space_ladder','active_space_ladder_boys','active_space_ladder_boys_large']:
        for p in sorted(Path('results/certificate_scaling',root).glob('h*/fixture.json')):
            f=json.loads(p.read_text());h=decode(f['hamiltonian'],f['modes'],4);_,_,s=partition(h,f['modes'],'quadratic',True)
            tails=[]
            for b in [F(0),F(1,625),F(1,10),F(1)]:
                kept,used=trim(h,b);tails.append({'coefficient_l1_budget':str(b),'used_bound':str(used),'retained_terms':len(kept),**graph(kept,f['modes'])})
            rows.append({'fixture':str(p),'modes':f['modes'],'natoms':f['natoms'],'hamiltonian_terms':len(h),
                'full_quadratic_matrix_entries':sum(n*n for n in s['original_dimensions']),
                'quotient_quadratic_matrix_entries':sum(n*n for n in s['symmetry_dimensions']),
                'largest_quotient_block':max(s['symmetry_dimensions']),'parity_generators':len(s['parity_masks']),
                'tails':tails})
    p=Path('results/certificate_scaling/structural_graph_diagnostics.json');p.write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps([{'fixture':r['fixture'],'entries':(r['full_quadratic_matrix_entries'],r['quotient_quadratic_matrix_entries']),
                      'graph_at_accuracy_budget':r['tails'][1]} for r in rows],indent=2))
if __name__=='__main__':main()
