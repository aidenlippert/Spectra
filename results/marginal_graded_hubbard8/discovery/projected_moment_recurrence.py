"""Exact four-site check of the nonorthogonal projected-moment recurrence.

This validates the algebra only. Both comparison routes still use sparse CAR
state actions; it is not a tensor-contraction implementation or scaling result.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.discovery.singlet_energy_replay import apply,gram
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_enlarged_schur import solve_positive
ROOT=Path(__file__).resolve().parents[1]
def mul(a,b):return [[sum(x*b[k][j] for k,x in enumerate(row)) for j in range(len(b[0]))] for row in a]
from experiments.marginal_symmetry_moments import projected_moments as recurrence
def main():
    h=json.loads((ROOT/'hamiltonian.json').read_text())
    h={'modes':8,'particles':4,'hamiltonian':[term for term in h['hamiltonian'] if all(mode<8 for _,mode in term['word'])]}
    o=DeterminantOracle(h);p=[sum(1<<(2*i+(i not in up)) for i in range(4)) for up in combinations(range(4),2)]
    V=[]
    for pairs in [((0,1),(2,3)),((0,3),(1,2))]:
        order=[s for pair in pairs for s in pair];phase=(-1)**sum(order[i]>order[j] for i in range(4) for j in range(i+1,4))
        v={}
        for s in p:
            spins=[(s>>(2*i))&1 for i in range(4)]
            if all(spins[i]!=spins[j] for i,j in pairs):v[s]=phase*(-1)**sum(not spins[i] for i,j in pairs)
        V.append(v)
    powers=V;M=[]
    for k in range(15):
        M.append(gram(V,powers))
        if k<14:powers=[apply(o,v) for v in powers]
    predicted=recurrence(M);W=[apply(o,v,set(p)) for v in V];powers=W;actual=[]
    for k in range(13):
        actual.append(gram(W,powers))
        if k<12:powers=[apply(o,v,set(p)) for v in powers]
    if predicted!=actual:raise ValueError('Projected moment recurrence mismatch')
    out=ROOT/'projected_moment_recurrence';out.mkdir(exist_ok=True)
    receipt={'accepted':True,'sites':4,'retained_dimension':2,'full_moment_orders':[0,14],
             'projected_moment_orders':[0,12],'all_entries_equal_exactly':True,
             'unique_action_states':len(o.cache),'referenced_determinants':o.referenced_state_count(),
             'scope':'Algebra validation only: sparse-action full moments versus direct Q-inserted moments. No tensor contraction or scalability claim.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
