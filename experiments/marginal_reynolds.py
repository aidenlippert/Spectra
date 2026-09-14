"""Full symmetry averaging of quartic SOS coefficient maps for matched models.

Discovery uses a small number of orbit representatives. Export expands their
averaged positive Gram matrices back to ordinary rational CAR certificates.
The floating symmetry reduction is never trusted for certificate acceptance.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time

import cvxpy as cp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

from experiments.marginal_coefficient import gram_map, export, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval
from experiments.marginal_dual_exact import component_charge, full_word_blocks
from experiments.marginal_orbit_quartic import symmetry_group, word_action
from experiments.marginal_stabilizer import split_blocks
from experiments.marginal_symbolic import add, canonical, multiplier_basis, number_shift, product, scale, transform


def invariant_rows(modes, degree):
    rows=[]
    for k in range(degree+1):
        sets=list(combinations(range(modes),k))
        for i,left in enumerate(sets):
            for right in sets[i:]:
                w=tuple((1,p) for p in left)+tuple((0,p) for p in right)
                if not any(component_charge(w,modes)):rows.append(w)
    return rows


def reynolds_matrix(rows, modes):
    """Full group average via signed orbits of generators, not enumeration."""
    lookup={w:i for i,w in enumerate(rows)};rr=[];cc=[];vv=[]
    half=modes//2;generators=[]
    for j in range(half-1):
        mapping=list(range(modes))
        for side in (0,half):mapping[side+j],mapping[side+j+1]=mapping[side+j+1],mapping[side+j]
        generators.append(mapping)
    generators.append([(i+half)%modes for i in range(modes)])
    unseen=set(range(len(rows)));chosen=[]
    while unseen:
        seed=min(unseen);signs={seed:1};queue=[seed];conflict=False
        for j in queue:
            for mapping in generators:
                image,sign=word_action(rows[j],mapping)
                if image not in lookup:
                    image=tuple((1,i) for c,i in image if not c)+tuple((0,i) for c,i in image if c)
                index=lookup[image];transport=int(sign)*signs[j]
                if index in signs:
                    if signs[index]!=transport:conflict=True
                else:signs[index]=transport;queue.append(index)
        unseen.difference_update(signs)
        # A stabilizer with sign -1 kills this invariant coefficient orbit.
        if conflict:continue
        chosen.append(seed)
        for i,si in signs.items():
            for j,sj in signs.items():rr.append(i);cc.append(j);vv.append(si*sj/len(signs))
    projection=sp.csr_matrix((vv,(rr,cc)),shape=(len(rows),len(rows)))
    return projection,chosen


def run(modes=10,degree=4,residual_penalty=False):
    started=time.monotonic();h=hopping_polynomial(modes,F(1,5));particles=modes//2
    rows=invariant_rows(modes,degree);lookup={w:i for i,w in enumerate(rows)}
    projection,chosen=reynolds_matrix(rows,modes);reduced=projection[chosen,:]
    full=full_word_blocks(modes,degree);split,_=split_blocks(modes,degree)
    representatives={}
    for i,block in enumerate(full):
        key=tuple(sorted(component_charge(block['words'][0],modes)))
        representatives.setdefault(key,i)
    indices=set(representatives.values())
    selected=[block for block in split if int(block['name'].split(':')[0]) in indices]
    basis=[p for p in multiplier_basis(modes,degree-1) if all(not any(component_charge(w,modes)) for w in p)]
    columns=[product(number_shift(modes,particles),p) for p in basis]
    a=sp.csr_matrix([[float(col.get(w,0)) for col in columns] for w in rows])
    pa=(reduced@a).toarray()
    _,r,piv=la.qr(pa,mode='economic',pivoting=True)
    rank=int(np.sum(abs(np.diag(r))>1e-10));piv=piv[:rank]
    # Full-group averaging of X is applied explicitly at export.
    x=cp.Variable(rank);b=cp.Variable();unit=np.array([float(w==()) for w in rows])
    variables=[];maps=[];dimensions=[]
    for block in selected:
        tr=block['basis_transform'];d=tr.shape[1]
        g=reduced@gram_map(block['words'],lookup)
        g=g@sp.kron(sp.csr_matrix(tr),sp.csr_matrix(tr),format='csr')
        maps.append(g);variables.append(cp.Variable((d,d),PSD=True));dimensions.append(d)
    expression=b*(reduced@unit)+pa[:,piv]@x+sp.hstack(maps,format='csr')@cp.hstack([
        cp.reshape(q,(d*d,),order='C') for q,d in zip(variables,dimensions)])
    rhs=np.array([float(h.get(w,0)) for w in rows]);objective=b
    if residual_penalty:
        residual=cp.Variable(len(chosen));weights=[]
        for i in chosen:
            # Reconstruct the full Hermitian polynomial's coefficient l1.
            row=projection.getrow(i);members=row.indices
            weights.append(sum(1 if tuple(j for c,j in rows[k] if c)==tuple(j for c,j in rows[k] if not c) else 2 for k in members))
        expression+=residual;objective=b-cp.norm1(cp.multiply(weights,residual))
    problem=cp.Problem(cp.Maximize(objective),[expression==reduced@rhs])
    build=time.monotonic()-started
    structure={'full_charge_blocks':len(full),'charge_orbits':len(indices),'psd_blocks':len(selected),
               'largest_block':max(dimensions),'psd_scalars':sum(d*(d+1)//2 for d in dimensions),
               'coefficient_rows_before':len(rows),'coefficient_rows_after':len(chosen),'multiplier_rank':rank}
    out=Path(__file__).resolve().parents[1]/'results/marginal_reynolds'/f'm{modes}_d{degree}_{"penalty" if residual_penalty else "exact"}'
    out.mkdir(parents=True,exist_ok=True);(out/'structure.json').write_text(json.dumps(structure,indent=2)+'\n')
    print(json.dumps(structure),flush=True);started=time.monotonic()
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=200)
    elapsed=time.monotonic()-started
    if b.value is None or x.value is None or any(q.value is None for q in variables):raise RuntimeError(problem.status)
    # Clip only the small numerical matrices, then average actual positive
    # matrices. The exact exported factors still bear all numerical errors.
    representative_grams={i:np.zeros((len(full[i]['words']),)*2) for i in indices}
    minimum=0.
    for block,q in zip(selected,variables):
        values,vectors=np.linalg.eigh((q.value+q.value.T)/2);minimum=min(minimum,float(values[0]))
        positive=(vectors*np.maximum(values,0))@vectors.T;tr=block['basis_transform']
        representative_grams[int(block['name'].split(':')[0])]+=tr@positive@tr.T
    group=symmetry_group(modes);grams=[np.zeros((len(block['words']),)*2) for block in full]
    charge_index={component_charge(block['words'][0],modes):i for i,block in enumerate(full)}
    word_indices=[{w:i for i,w in enumerate(block['words'])} for block in full]
    for i,gram in representative_grams.items():
        for mapping in group:
            images=[word_action(w,mapping) for w in full[i]['words']]
            target=charge_index[component_charge(images[0][0],modes)]
            permutation=[word_indices[target][w] for w,s in images];signs=np.array([s for w,s in images],dtype=float)
            grams[target][np.ix_(permutation,permutation)]+=gram*signs[:,None]*signs[None,:]/len(group)
    raw_x=add(*(scale(basis[int(j)],float(value)) for j,value in zip(piv,x.value)))
    averaged_x={}
    for mapping in group:
        for w,c in transform(raw_x,mapping).items():averaged_x[w]=averaged_x.get(w,0)+c/len(group)
    # Export accepts a polynomial-valued basis; one basis element suffices.
    solution={'b':float(b.value),'grams':grams,'x':np.array([1.]),'basis':[canonical(averaged_x)],
              'status':problem.status,'coefficient_rows':len(chosen),'gram_dimensions':dimensions,
              'coefficient_map_nonzeros':sum(g.nnz for g in maps),'multiplier_basis_dimension':rank,
              'build_seconds':build,'solve_seconds':elapsed}
    # Make X rational before export; no floating polynomial coefficients may
    # enter the exact verifier's algebra.
    solution['basis']=[{w:F(round(c*10**12),10**12) for w,c in averaged_x.items() if round(c*10**12)}]
    certificate,receipt=export(h,modes,particles,full,solution,denominator=10**9,operator_degree=degree)
    certificate['variational_upper']=symmetric_upper(modes,F(1,5));receipt.update(verify_interval(certificate))
    receipt.update(structure);receipt.update({'numerical_objective':float(problem.value),'minimum_raw_gram_eigenvalue':minimum,'residual_penalty':residual_penalty})
    (out/'certificate.json').write_text(json.dumps(certificate)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('lower_float','width_float','numerical_objective','minimum_raw_gram_eigenvalue','solve_seconds')}),flush=True)
    return certificate,receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--modes',type=int,default=10);parser.add_argument('--degree',type=int,default=4);parser.add_argument('--penalty',action='store_true')
    args=parser.parse_args();run(args.modes,args.degree,args.penalty)
