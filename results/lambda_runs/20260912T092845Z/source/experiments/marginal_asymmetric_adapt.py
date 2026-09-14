"""Adapt cubic squares around a compact quartic orbit certificate.

Only global left/right exchange and individual pair charges are used for the
new Hamiltonian. Quartic orbit directions remain positive without requiring
flavor permutation symmetry of H. No many-body vector enters the lower solve.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import time

import cvxpy as cp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

from experiments.marginal_coefficient import gram_map
from experiments.marginal_dual_exact import full_word_blocks,component_charge
from experiments.marginal_reynolds import invariant_rows,reynolds_matrix
from experiments.marginal_orbit_quartic import word_action,symmetry_group
from experiments.marginal_orbit_certificate import average_canonical_polynomial,normalize_scalar_residual,verify
from experiments.marginal_symbolic import adj,add,decode,encode,multiplier_basis,number_shift,product,scale,transform
from experiments.marginal_symmetry_transfer import hopping_perturbation,integer_upper,replay

ROOT=Path(__file__).resolve().parents[1]


def exchange_subspaces(words,modes):
    mapping=[(i+modes//2)%modes for i in range(modes)]
    lookup={w:i for i,w in enumerate(words)};seen=set();columns={1:[],-1:[]}
    for j,w in enumerate(words):
        if j in seen:continue
        image,sign=word_action(w,mapping);k=lookup[image];seen.update((j,k))
        if k==j:
            v=np.zeros(len(words));v[j]=1;columns[int(sign)].append(v)
        else:
            for parity in (1,-1):
                v=np.zeros(len(words));v[j]=1/np.sqrt(2);v[k]=parity*sign/np.sqrt(2)
                columns[parity].append(v)
    return [np.column_stack(v) for v in columns.values() if v]


def exchange_projection(rows,modes):
    mapping=[(i+modes//2)%modes for i in range(modes)];lookup={w:i for i,w in enumerate(rows)}
    rr=[];cc=[];vv=[];weights=[];seen=set()
    for j,w in enumerate(rows):
        if j in seen:continue
        image,sign=word_action(w,mapping)
        if image not in lookup:image=tuple((1,i) for c,i in image if not c)+tuple((0,i) for c,i in image if c)
        k=lookup[image];seen.update((j,k))
        if k==j and sign==-1:continue
        row=len(weights);members=[j] if k==j else [j,k]
        for i in members:rr.append(row);cc.append(i);vv.append(float((1 if i==j else sign)/len(members)))
        weights.append(sum(1 if tuple(i for c,i in rows[a] if c)==tuple(i for c,i in rows[a] if not c) else 2 for a in members))
    return sp.csr_matrix((vv,(rr,cc)),shape=(len(weights),len(rows))),np.array(weights)


def prepare(source,aggregate=False):
    started=time.monotonic();base=json.loads(Path(source).read_text());modes=base['modes'];particles=base['particles']
    if modes%2 or modes<4:raise ValueError('Expected even matched mode count')
    h=decode(base['hamiltonian'],modes,4);mapping=[(i+modes//2)%modes for i in range(modes)]
    if transform(h,mapping)!=h or any(any(component_charge(w,modes)) for w in h):raise ValueError('Source breaks pair charges or exchange')
    rows=invariant_rows(modes,3 if aggregate else 4);lookup={w:i for i,w in enumerate(rows)};projection,weights=exchange_projection(rows,modes)
    basis=[p for p in multiplier_basis(modes,2 if aggregate else 3) if all(not any(component_charge(w,modes)) for w in p)]
    ideal=[product(number_shift(modes,particles),p) for p in basis]
    a=projection@np.array([[float(col.get(w,0)) for col in ideal] for w in rows])
    _,r,piv=la.qr(a,mode='economic',pivoting=True);rank=int(np.sum(abs(np.diag(r))>1e-10));piv=piv[:rank]
    blocks=[];maps=[]
    for block in full_word_blocks(modes,3):
        g=projection@gram_map(block['words'],lookup)
        for tr in exchange_subspaces(block['words'],modes):
            blocks.append({'words':block['words'],'transform':tr});st=sp.csr_matrix(tr)
            maps.append(g@sp.kron(st,st,format='csr'))
    seeds=[];columns=[]
    full_average=None
    if not aggregate and {tuple(g) for g in base['permutations']}==set(symmetry_group(modes)):
        full_average=reynolds_matrix(rows,modes)[0]
    for item in ([] if aggregate else base['orbit_squares']):
        p=decode(item['polynomial'],modes,4);norm=max(abs(c) for c in p.values());p=scale(p,1/norm)
        square=product(adj(p),p)
        if full_average is None:square=average_canonical_polynomial(square,base['permutations'])
        column=np.array([float(square.get(w,0)) for w in rows])
        if full_average is not None:column=full_average@column
        seeds.append(p);columns.append(projection@column)
    source_receipt=verify(base) if aggregate else None
    if aggregate:
        atom=dict(h);atom[()]=atom.get((),F(0))-F(source_receipt['lower'])
        columns=[projection@np.array([float(atom.get(w,0)) for w in rows])]
    return {'base':base,'modes':modes,'particles':particles,'h':h,'rows':rows,'projection':projection,'weights':weights,
            'basis':[basis[int(j)] for j in piv],'a':a[:,piv],'blocks':blocks,'maps':maps,'seeds':seeds,'columns':np.column_stack(columns),
            'unit':projection@np.array([float(w==()) for w in rows]),'build_seconds':time.monotonic()-started,
            'aggregate':aggregate,'source_receipt':source_receipt}


def run(epsilon=F(1,1000),source=None,max_seconds=120,penalty=False,aggregate=False):
    source=ROOT/'results/marginal_distill_direct/compact_certificate.json' if source is None else Path(source)
    data=prepare(source,aggregate);modes=data['modes'];particles=data['particles'];delta,_,_=hopping_perturbation(modes,epsilon);h=add(data['h'],delta)
    dims=[block['transform'].shape[1] for block in data['blocks']];n=data['projection'].shape[0]
    b=cp.Variable();x=cp.Variable(len(data['basis']));alpha=cp.Variable(data['columns'].shape[1],nonneg=True);residual=cp.Variable(n) if penalty else np.zeros(n)
    grams=[cp.Variable((d,d),PSD=True) for d in dims]
    expression=b*data['unit']+data['a']@x+data['columns']@alpha+residual
    expression+=sp.hstack(data['maps'],format='csr')@cp.hstack([cp.reshape(q,(d*d,),order='C') for q,d in zip(grams,dims)])
    equality=expression==data['projection']@np.array([float(h.get(w,0)) for w in data['rows']])
    constraints=[equality,alpha<=100]
    if penalty:constraints.append(residual[0]==0)
    problem=cp.Problem(cp.Maximize(b-cp.norm1(cp.multiply(data['weights'],residual)) if penalty else b),constraints)
    out=ROOT/'results/marginal_asymmetric_adapt'/(str(epsilon).replace('/','_')+('_atom' if aggregate else '')+('_penalty' if penalty else '_exact'));out.mkdir(parents=True,exist_ok=True)
    structure={'epsilon':str(epsilon),'source':str(source),'coefficient_rows':n,'full_rows':len(data['rows']),'multiplier_rank':len(data['basis']),
               'cubic_blocks':len(dims),'largest_cubic_block':max(dims),'cubic_psd_scalars':sum(d*(d+1)//2 for d in dims),
               'quartic_scalar_directions':data['columns'].shape[1],'build_seconds':data['build_seconds'],'solver_time_limit':max_seconds,'residual_penalty':penalty,
               'certified_operator_atom':aggregate,'maximum_scalar_weight':100}
    (out/'structure.json').write_text(json.dumps(structure,indent=2)+'\n');print(json.dumps(structure),flush=True)
    started=time.monotonic()
    try:
        problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=150,time_limit=max_seconds)
    except cp.error.SolverError as error:
        status={'status':'solver_error','solve_seconds':time.monotonic()-started,'error':str(error)}
        (out/'solver_status.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status),flush=True);return status
    elapsed=time.monotonic()-started
    status={'status':problem.status,'solve_seconds':elapsed,'numerical_objective':float(problem.value) if problem.value is not None else None}
    (out/'solver_status.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status),flush=True)
    if b.value is None or x.value is None or alpha.value is None or any(q.value is None for q in grams):return status
    def rational(value):return F(round(float(value)*10**12),10**12)
    direct=[];minimum_eigenvalue=0.
    for block,q in zip(data['blocks'],grams):
        values,vectors=np.linalg.eigh((q.value+q.value.T)/2)
        minimum_eigenvalue=min(minimum_eigenvalue,float(values[0]))
        for j,value in enumerate(values):
            if value<=0:continue
            v=block['transform']@vectors[:,j]*np.sqrt(value)
            p={w:rational(c) for w,c in zip(block['words'],v) if rational(c)}
            if p:direct.append({'polynomial':encode(p),'weight':'1'})
    orbit=[{'polynomial':encode(p),'weight':str(rational(max(0,w)))} for p,w in zip(data['seeds'],alpha.value) if rational(max(0,w))]
    multiplier=add(*(scale(p,rational(c)) for p,c in zip(data['basis'],x.value)))
    multiplier=average_canonical_polynomial(multiplier,[list(range(modes)),[(i+modes//2)%modes for i in range(modes)]])
    exported_b=rational(b.value)
    if aggregate:
        weight=rational(max(0,alpha.value[0]));base=data['base']
        orbit=[dict(item,weight=str(weight*F(item['weight']))) for item in base['orbit_squares']]
        direct.extend(dict(item,weight=str(weight*F(item['weight']))) for item in base.get('direct_squares',[]))
        multiplier=add(multiplier,scale(decode(base['number_multiplier'],modes,6),weight))
        exported_b+=weight*(F(base['b'])-F(data['source_receipt']['lower']))
    certificate={'modes':modes,'particles':particles,'operator_degree':4,'hamiltonian':encode(h),'number_multiplier':encode(multiplier),
                 'b':str(exported_b),'permutations':data['base']['permutations'],'orbit_squares':orbit,'direct_squares':direct}
    certificate,_=normalize_scalar_residual(certificate)
    certificate['independent_upper']=integer_upper(h,modes,particles);receipt=replay(certificate)
    receipt['minimum_raw_gram_eigenvalue']=minimum_eigenvalue
    receipt.update(structure);receipt.update(status)
    (out/'certificate.json').write_text(json.dumps(certificate)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'lower':receipt['lower_float'],'upper':receipt['upper_float'],'width':receipt['width_float']}),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--epsilon',default='1/1000');parser.add_argument('--source',type=Path);parser.add_argument('--time-limit',type=float,default=120)
    parser.add_argument('--penalty',action='store_true')
    parser.add_argument('--aggregate',action='store_true')
    args=parser.parse_args();run(F(args.epsilon),args.source,args.time_limit,args.penalty,args.aggregate)
