"""Distill an accepted symmetric certificate into a few orbits of squares.

Directions come from an existing factor or orbit certificate. A linear program
chooses weights; an independent rational checker accepts the new certificate.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import linprog
import scipy.linalg as la
import scipy.sparse as sp

from experiments.marginal_dual_exact import component_charge
from experiments.marginal_reynolds import invariant_rows, reynolds_matrix
from experiments.marginal_orbit_quartic import symmetry_group
from experiments.marginal_symbolic import adj, add, decode, encode, multiplier_basis, number_shift, product, scale, transform


def accept_compact(certificate):
    from experiments.marginal_orbit_certificate import verify
    from experiments.marginal_collective import hopping_polynomial, upper_bound
    for item in certificate['orbit_squares']:
        p=decode(item['polynomial'],10,4)
        item['polynomial']=encode({w:F(round(a*10**14),10**14) for w,a in p.items() if round(a*10**14)})
    receipt=verify(certificate)
    if certificate['modes']!=10 or certificate['particles']!=5 or decode(certificate['hamiltonian'],10,4)!=hopping_polynomial(10,F(1,5)):
        raise ValueError('Upper witness requires the matched ten-mode Hamiltonian')
    upper=upper_bound(10,F(1,5),certificate['variational_upper']['amplitudes'])['upper']
    receipt.update({'upper':upper,'width':str(F(upper)-F(receipt['lower']))})
    receipt['width_float']=float(F(receipt['width']))
    if F(receipt['width'])<0:raise ValueError('Inconsistent interval')
    return certificate,receipt


def candidate_seeds(original):
    modes=original['modes'];degree=original['operator_degree'];polynomials=[]
    if 'orbit_squares' in original:
        if {tuple(g) for g in original['permutations']}!=set(symmetry_group(modes)):
            raise ValueError('Distillation expects the full matched symmetry group')
        polynomials=[decode(item['polynomial'],modes,degree) for item in original['orbit_squares'] if F(item['weight'])>0]
    else:
        representatives={};denominator=original['denominator']
        for block in original['blocks']:
            words=[tuple(map(tuple,w)) for w in block['words']]
            key=tuple(sorted(component_charge(words[0],modes)))
            representatives.setdefault(key,(words,block['factor']))
        for words,factors in representatives.values():
            polynomials.extend({w:F(c,denominator) for w,c in zip(words,factor) if c} for factor in factors)
    seeds=[]
    for p in polynomials:
        if sum(float(c)**2 for c in p.values())<1e-14:continue
        biggest=max(abs(c) for c in p.values())
        seeds.append({w:c/biggest for w,c in p.items()})
    return seeds


def run(source=None,out_name='marginal_distill'):
    root=Path(__file__).resolve().parents[1]
    source=root/'results/marginal_reynolds/m10_d4_exact/certificate.json' if source is None else Path(source)
    original=json.loads(source.read_text())
    modes=original['modes'];particles=original['particles'];degree=original['operator_degree']
    h=decode(original['hamiltonian'],modes,4);rows=invariant_rows(modes,degree)
    projection,chosen=reynolds_matrix(rows,modes);reduced=projection[chosen,:]
    seeds=candidate_seeds(original);columns=[]
    for p in seeds:
        square=product(adj(p),p)
        columns.append(reduced@np.array([float(square.get(w,0)) for w in rows]))
    basis=[p for p in multiplier_basis(modes,degree-1) if all(not any(component_charge(w,modes)) for w in p)]
    ideal=[product(number_shift(modes,particles),p) for p in basis]
    a=reduced@np.array([[float(p.get(w,0)) for p in ideal] for w in rows])
    _,r,piv=la.qr(a,mode='economic',pivoting=True);rank=int(np.sum(abs(np.diag(r))>1e-10));piv=piv[:rank]
    unit=reduced@np.array([float(w==()) for w in rows])
    matrix=np.column_stack([unit,a[:,piv],*columns]);rhs=reduced@np.array([float(h.get(w,0)) for w in rows])
    objective=np.zeros(matrix.shape[1]);objective[0]=-1
    started=time.monotonic()
    result=linprog(objective,A_eq=sp.csr_matrix(matrix),b_eq=rhs,
                   bounds=[(None,None)]*(1+rank)+[(0,100)]*len(seeds),method='highs',
                   options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not result.success:raise RuntimeError(result.message)
    group=symmetry_group(modes)
    raw_x=add(*(scale(basis[int(j)],F(round(float(c)*10**14),10**14)) for j,c in zip(piv,result.x[1:1+rank])))
    averaged_x=scale(add(*(transform(raw_x,mapping) for mapping in group)),F(1,len(group)))
    squares=[{'polynomial':encode(p),'weight':str(F(round(float(weight)*10**14),10**14))}
             for p,weight in zip(seeds,result.x[1+rank:]) if weight>1e-10]
    certificate={'modes':modes,'particles':particles,'operator_degree':degree,'hamiltonian':encode(h),
                 'number_multiplier':encode(averaged_x),'b':str(F(round(float(result.x[0])*10**14),10**14)),
                 'permutations':[list(g) for g in group],'orbit_squares':squares,
                 'variational_upper':original['variational_upper']}
    out=root/'results'/out_name/'bounded';out.mkdir(parents=True,exist_ok=True)
    (out/'certificate.json').write_text(json.dumps(certificate)+'\n')
    proposal={'source':str(source),'candidate_rays':len(seeds),'selected_rays':len(squares),'rows':len(chosen),'multiplier_rank':rank,
              'numerical_b':float(result.x[0]),'lp_seconds':time.monotonic()-started,
              'numerical_matching_max_error':float(np.max(abs(matrix@result.x-rhs))),
              'seed_monomials':sum(len(s['polynomial']) for s in squares)}
    (out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n');print(json.dumps(proposal),flush=True)
    compact,receipt=accept_compact(certificate)
    (out.parent/'compact_certificate.json').write_text(json.dumps(compact)+'\n')
    (out.parent/'compact_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'accepted_lower':receipt['lower_float'],'accepted_width':receipt['width_float']}),flush=True)
    return certificate,proposal


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path);parser.add_argument('--output-name',default='marginal_distill')
    args=parser.parse_args();run(args.source,args.output_name)
