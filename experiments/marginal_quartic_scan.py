"""Coefficient-space scan of quartic charged operator sectors."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

from experiments.marginal_coefficient import dagger, export, gram_map, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval
from experiments.marginal_compression import charge, hopping_components
from experiments.marginal_cubic_probe import complete_cubic
from experiments.marginal_adaptive import assemble_and_solve
from experiments.marginal_symbolic import add, adj, mono, multiplier_basis, number_shift, product, scale


def pair_polynomials(modes):
    half=modes//2
    hops=[mono(((1,i),(0,i+half))) for i in range(half)]
    return [product(hops[i],hops[j]) for i,j in combinations(range(half),2)]


def probe_block(modes,kind):
    forward=pair_polynomials(modes)
    if kind=='antisymmetric':
        polynomials=[add(a,scale(adj(a),-1)) for a in forward]
    elif kind=='span':
        half=modes//2
        polynomials=[mono(())]
        polynomials += [mono(((1,i),(0,j))) for pair in ((k,k+half) for k in range(half)) for i in pair for j in pair]
        hops=[mono(((1,i),(0,i+half))) for i in range(half)]
        polynomials += [product(a,b) for i,j in combinations(range(half),2)
                        for a in (hops[i],adj(hops[i])) for b in (hops[j],adj(hops[j]))]
    elif kind == 'density_dressed':
        half = modes // 2
        hops = [mono(((1,i),(0,i+half))) for i in range(half)]
        densities = [mono(((1,i),(0,i))) for i in range(modes)]
        polynomials = [product(n, h) for n in densities for h in hops]
        polynomials += [adj(p) for p in polynomials]
    elif kind == 'quartet':
        polynomials = [mono(tuple((0,i) for i in q)) for q in combinations(range(modes), 4)]
        polynomials += [adj(p) for p in polynomials]
    elif kind == 'charged2':
        # 1-create/3-annihilator quartics, with their adjoints. Add quadratic
        # pair annihilators in the same charge sectors after splitting.
        polynomials = []
        for c in range(modes):
            for a in combinations([i for i in range(modes) if i != c], 3):
                polynomials.append(mono(((1,c),)+tuple((0,i) for i in a)))
        polynomials += [adj(p) for p in polynomials]
        polynomials += [mono(tuple((0,i) for i in a)) for a in combinations(range(modes),2)]
        polynomials += [adj(p) for p in polynomials[-len(list(combinations(range(modes),2))):]]
    elif kind == 'neutral4':
        polynomials=[]
        for c in combinations(range(modes),2):
            for a in combinations([i for i in range(modes) if i not in c],2):
                polynomials.append(mono(tuple((1,i) for i in c)+tuple((0,i) for i in a)))
    else:
        raise ValueError('Unknown pair probe')
    words=sorted({tuple((int(c),i) for c,i in w) for p in polynomials for w in p})
    transform=np.array([[float(p.get(w,0)) for p in polynomials] for w in words])
    return {'name':kind,'words':words},transform


def solve(modes=10,t=F(1,5),kind='antisymmetric',solver='SCS'):
    start=time.monotonic()
    h=hopping_polynomial(modes,t); groups=hopping_components(h,modes)
    membership={i:g for g,group in enumerate(groups) for i in group}
    invariant=lambda w:not any(charge(w,membership,len(groups)))
    rows=[]
    for k in range(5):
        sets=list(combinations(range(modes),k))
        for i,left in enumerate(sets):
            for right in sets[i:]:
                w=tuple((1,p) for p in left)+tuple((0,p) for p in right)
                if invariant(w): rows.append(w)
    lookup={w:i for i,w in enumerate(rows)}
    basis=[p for p in multiplier_basis(modes,3) if all(invariant(w) for w in p)]
    shift=number_shift(modes,modes//2)
    ri=[]; ci=[]; values=[]
    for j,p in enumerate(basis):
        for w,c in product(shift,p).items():
            if w in lookup:ri.append(lookup[w]);ci.append(j);values.append(float(c))
    a=sp.csr_matrix((values,(ri,ci)),shape=(len(rows),len(basis)))
    b=cp.Variable(); x=cp.Variable(len(basis))
    unit=np.zeros(len(rows));unit[lookup[()]]=1
    expression=b*unit+a@x
    blocks=complete_cubic(h,modes)
    transforms=[None]*len(blocks)
    if kind!='control':
        block,transform=probe_block(modes,kind)
        blocks.append(block);transforms.append(transform)
    variables=[]; dimensions=[]; nnz=0
    for block,transform in zip(blocks,transforms):
        words=block['words']; dimension=len(words) if transform is None else transform.shape[1]
        q=cp.Variable((dimension,dimension),PSD=True)
        expanded=q if transform is None else transform@q@transform.T
        g=gram_map(words,lookup)
        expression += g@cp.reshape(expanded,(len(words)**2,),order='C')
        variables.append(q);dimensions.append(dimension);nnz+=g.nnz
    problem=cp.Problem(cp.Maximize(b),[expression==np.array([float(h.get(w,0)) for w in rows])])
    build=time.monotonic()-start;start=time.monotonic()
    if solver=='SCS':
        problem.solve(solver='SCS',eps=1e-9,max_iters=100000)
    elif solver=='CLARABEL':
        problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=200)
    else:
        raise ValueError('Unsupported probe solver')
    elapsed=time.monotonic()-start
    if b.value is None or x.value is None or any(q.value is None for q in variables):
        raise RuntimeError(f'No proposal: {problem.status}')
    grams=[q.value if tr is None else tr@q.value@tr.T for q,tr in zip(variables,transforms)]
    solution={'b':float(b.value),'x':np.asarray(x.value),'basis':basis,'grams':grams,
              'status':problem.status,'coefficient_rows':len(rows),'gram_dimensions':dimensions,
              'coefficient_map_nonzeros':nnz,'multiplier_basis_dimension':len(basis),
              'build_seconds':build,'solve_seconds':elapsed}
    cert,receipt=export(h,modes,modes//2,blocks,solution,operator_degree=4)
    cert['variational_upper']=symmetric_upper(modes,t)
    receipt.update(verify_interval(cert));receipt['probe']=kind
    receipt['solver']=solver
    receipt['solver_eps']=1e-9;receipt['max_iters']=100000 if solver=='SCS' else 200
    receipt['iterations']=problem.solver_stats.num_iters
    return cert,receipt

def solve_adaptive(modes=10, t=F(1,5), kind='quartet'):
    h=hopping_polynomial(modes,t); groups=hopping_components(h,modes)
    member={i:g for g,grp in enumerate(groups) for i in grp}
    def q(w):
        z=[0]*len(groups)
        for c,i in w:z[member[i]] += 2*c-1
        return tuple(z)
    base=complete_cubic(h,modes)
    extra,_=probe_block(modes,kind) if kind != 'union' else ({'words':[]},None)
    from collections import defaultdict
    sectors=defaultdict(list)
    for w in extra['words']: sectors[q(w)].append(w)
    if kind == 'union':
        extras=[probe_block(modes,k)[0] for k in ('charged2','quartet','neutral4')]
        sectors=defaultdict(list)
        for ex in extras:
            for w in ex['words']: sectors[q(w)].append(w)
    blocks=base+[{'name':f'{kind}:{charge}', 'words':list(dict.fromkeys(ws))} for charge,ws in sectors.items()]
    proposal,data=assemble_and_solve(h,modes,modes//2,blocks,degree=4,groups=groups)
    cert,receipt=export(h,modes,modes//2,blocks,proposal,operator_degree=4)
    cert['variational_upper']=symmetric_upper(modes,t); receipt.update(verify_interval(cert))
    receipt.update({'probe':kind,'solver':'CLARABEL','quartic_sectors':list(sectors),'adaptive':True})
    return cert,receipt

def solve_complete():
    from experiments.marginal_dual_exact import full_word_blocks
    h=hopping_polynomial(10,F(1,5)); groups=hopping_components(h,10)
    blocks=full_word_blocks(10,4)
    proposal,data=assemble_and_solve(h,10,5,blocks,degree=4,groups=groups)
    cert,receipt=export(h,10,5,blocks,proposal,operator_degree=4)
    cert['variational_upper']=symmetric_upper(10,F(1,5)); receipt.update(verify_interval(cert))
    receipt.update({'probe':'complete4','solver':'CLARABEL','block_count':len(blocks)})
    out=Path(__file__).resolve().parents[1]/'results/marginal_quartic_scan';out.mkdir(exist_ok=True)
    (out/'complete4.json').write_text(json.dumps(cert)+'\n')
    (out/'complete4_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',choices=('control','antisymmetric','span','density_dressed','quartet','charged2','neutral4','union','complete4'),default='charged2')
    parser.add_argument('--solver',choices=('SCS','CLARABEL'),default='SCS')
    args=parser.parse_args()
    if args.probe == 'complete4': solve_complete(); return
    cert,receipt=solve_adaptive(kind=args.probe)
    out=Path(__file__).resolve().parents[1]/'results/marginal_quartic_scan';out.mkdir(exist_ok=True)
    name=args.probe+('_clarabel' if args.solver=='CLARABEL' else '')
    (out/f'{name}.json').write_text(json.dumps(cert)+'\n')
    (out/f'{name}_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__': main()
