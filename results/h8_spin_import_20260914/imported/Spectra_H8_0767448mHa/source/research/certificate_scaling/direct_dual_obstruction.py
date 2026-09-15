"""Direct LP dual obstruction for the complete quadratic square cone (H4).

This is a family-specific diagnostic.  It does not use FCI or a sector matrix.
Numerical LP output is never called a proof: the rationalized dual is checked
against every exact ideal and square constraint before being reported.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.marginal_coefficient import dictionaries, dagger
from experiments.marginal_symbolic import add, canonical, multiplier_basis, mono, number_shift, word_product


def square(words, signs=(1, 1)):
    out = {}
    for i, left in enumerate(words):
        for j, right in enumerate(words):
            c = signs[i]*signs[j]
            for w, s in word_product(dagger(left), right): out[w] = out.get(w,F(0))+F(c*s)
    return {w:c for w,c in out.items() if c}


def addmap(dst, src, factor=F(1)):
    for w,c in src.items(): dst[w] = dst.get(w,F(0))+factor*c


def run():
    M,N=8,4
    blocks=dictionaries(M,'quadratic'); words=[]
    for b in blocks: words.extend(tuple(tuple(x) for x in w) for w in b['words'])
    atoms=[square([w]) for w in words]
    for i,j in combinations(range(len(words)),2):
        atoms.append(square([words[i],words[j]]))
        atoms.append(square([words[i],words[j]],(1,-1)))
    # Number-ideal constraints are (Nhat-N)*X with X body <=2.
    shift=number_shift(M,N); ideals=[]
    for x in multiplier_basis(M,2):
        p={}
        for a,ca in shift.items():
            for b,cb in x.items():
                for w,s in word_product(a,b): p[w]=p.get(w,F(0))+ca*cb*s
        ideals.append({w:c for w,c in p.items() if c})
    h={tuple(tuple(t) for t in q['word']):F(q['coefficient']) for q in json.loads((Path(__file__).resolve().parents[2]/'results/marginal_molecule_stress/h4_square_degree3_certificate.json').read_text())['hamiltonian']}
    universe={()}
    for p in atoms+ideals+[h]: universe.update(p)
    universe=sorted(universe,key=lambda w:(len(w),w)); idx={w:i for i,w in enumerate(universe)}; n=len(universe)
    def vec(p): return np.array([float(p.get(w,0)) for w in universe])
    Aeq=np.zeros((1+len(ideals),n)); Aeq[0,idx[()]]=1
    for k,p in enumerate(ideals,1): Aeq[k]=vec(p)
    Aub=np.zeros((len(atoms),n))
    for k,p in enumerate(atoms): Aub[k]=-vec(p)
    res=linprog(vec(h), A_ub=csr_matrix(Aub), b_ub=np.zeros(len(atoms)), A_eq=csr_matrix(Aeq), b_eq=np.r_[1.,np.zeros(len(ideals))], bounds=[(-1,1)]*n, method='highs')
    y=[F(float(v)).limit_denominator(10**6) for v in res.x] if res.success else []
    exact_ideal=max((abs(sum((c*y[idx[w]] for w,c in p.items()),F(0))) for p in ideals),default=F(0))
    exact_atom_min=min((sum((c*y[idx[w]] for w,c in p.items()),F(0)) for p in atoms),default=F(0))
    exact_const=y[idx.get((),0)] if y else None; norm=max((abs(v) for v in y),default=F(0))
    exact_objective=sum((c*y[idx[w]] for w,c in h.items()),F(0)) if y else None
    valid=bool(res.success and exact_const==1 and exact_ideal==0 and exact_atom_min>=0 and norm<=1)
    out={'scope':'complete quadratic singleton and signed-pair squares, M=8,N=4','words':len(words),'atoms':len(atoms),'universe_words':n,'lp_success':bool(res.success),'lp_message':res.message,'numeric_objective':float(res.fun) if res.success else None,'exact_objective':str(exact_objective),'rationalized_valid':valid,'exact_checks':{'constant':str(exact_const),'max_ideal_abs':str(exact_ideal),'min_atom':str(exact_atom_min),'norm_inf':str(norm)},'dual_words':[{"word":[list(x) for x in w],"value":str(y[k])} for k,w in enumerate(universe) if y and y[k]],'cost':{'lp_variables':n,'ineq':len(atoms),'equalities':1+len(ideals)}}
    dest=Path(__file__).resolve().parents[2]/'results/certificate_scaling/direct_dual_obstruction'; dest.mkdir(parents=True,exist_ok=True); (dest/'receipt.json').write_text(json.dumps(out,indent=2)+'\n'); return out

if __name__=='__main__': print(json.dumps(run(),indent=2))
