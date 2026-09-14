"""Bounded direct sparse CAR-SOS atom discovery (no full Gram matrices)."""
import argparse, json, time
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix
from fractions import Fraction as F
from experiments.marginal_coefficient import coefficient_rows, gram_map, hopping_model
from experiments.marginal_symbolic import multiplier_basis, number_shift, product

def words(m):
    return [((0,i),) for i in range(m)] + [((0,j),(0,i)) for i,j in combinations(range(m),2)]

def build_atoms(modes, max_atoms=256):
    ws=words(modes); atoms=[]
    # Rank-one squares of monomials and normalized +/- two-word directions.
    for w in ws: atoms.append((w, w, 1.0))
    charge=lambda w:sum(2*c-1 for c,_ in w)
    for i,j in combinations(range(len(ws)),2):
        if charge(ws[i]) != charge(ws[j]): continue
        if len(atoms)>=max_atoms: break
        atoms.extend([(ws[i],ws[j],1.0),(ws[i],ws[j],-1.0)])
        if len(atoms)>=max_atoms: break
    return atoms[:max_atoms]

def run(modes=4, max_atoms=128, residual_penalty=1e3):
    started=time.perf_counter(); h=hopping_model(modes,F(1,5),False)
    rows=coefficient_rows(modes); lookup={w:i for i,w in enumerate(rows)}; atoms=build_atoms(modes,max_atoms)
    maps=[]; kept=[]
    for a,b,s in atoms:
        ws=[a,b] if a!=b else [a]
        v=np.array([1.,s]) if a!=b else np.array([1.])
        g=gram_map(ws,lookup)
        col=np.asarray(g @ np.outer(v,v).reshape(-1,order='C')).ravel()
        if np.linalg.norm(col)>0: maps.append(col); kept.append((a,b,s))
    basis=multiplier_basis(modes); shift=number_shift(modes,modes//2); ideal=[]
    for q in basis:
        p=product(shift,q); ideal.append(np.array([float(p.get(w,0)) for w in rows]))
    unit=np.zeros(len(rows)); unit[lookup[()]]=1.0
    A=np.column_stack([unit]+ideal+maps); rhs=np.array([float(h.get(w,0)) for w in rows])
    nvar=A.shape[1]; nb=len(ideal); na=len(maps); slack=len(rows)
    # variables [b, x_free, weights>=0, residual+>=0, residual->=0]
    c=np.r_[ -1., np.zeros(nb), np.zeros(na), np.full(2*slack,residual_penalty)]
    eq=np.column_stack([A, np.eye(slack), -np.eye(slack)])
    bounds=[(None,None)]+[(None,None)]*nb+[(0,None)]*na+[(0,None)]*(2*slack)
    sol=linprog(c,A_eq=csc_matrix(eq),b_eq=rhs,bounds=bounds,method='highs')
    out=Path(__file__).resolve().parents[2]/'results/certificate_scaling/direct_discovery';out.mkdir(parents=True,exist_ok=True)
    rec={'scope':'direct rank-one CAR SOS atom LP; numerical proposal only','modes':modes,'rows':len(rows),'candidate_atoms':len(atoms),'nonzero_maps':na,'ideal_columns':nb,'max_atoms':max_atoms,'status':sol.message,'success':bool(sol.success),'wall_seconds':time.perf_counter()-started}
    if sol.success:
        represented=A@sol.x[:1+nb+na]; residual=represented-rhs; weights=sol.x[1+nb:1+nb+na]; rec.update({'b':float(sol.x[0]),'residual_l1':float(np.sum(abs(residual))),'residual_max':float(np.max(abs(residual))),'active_atoms':int(np.count_nonzero(weights>1e-8)),'active_weight_sum':float(weights.sum()),'atom_support':[kept[i][0] for i,w in enumerate(weights) if w>1e-8]})
    (out/f'm{modes}_a{max_atoms}.json').write_text(json.dumps(rec,indent=2)+'\n'); return rec

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--modes',type=int,default=4);p.add_argument('--max-atoms',type=int,default=128);a=p.parse_args();print(json.dumps(run(a.modes,a.max_atoms),indent=2))
