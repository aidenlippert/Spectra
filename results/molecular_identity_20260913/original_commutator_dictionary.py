"""H-only linear/commutator polynomial dictionaries with exact CAR export.

The full quadratic baseline is retained. Cubic parts of [H,a_i] enrich its
linear charge blocks. Every coefficient of every square is matched, including
degree six. Small Gram dimension does not hide polynomial-product work.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import argparse,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
import cvxpy as cp
from experiments.marginal_symbolic import canonical,hermitian,decode,encode,add,scale,mono,product,number_shift,multiplier_basis,verify
from experiments.marginal_hunt_car import adj
from research.certificate_scaling.adaptive_block_discovery import partition
from research.certificate_scaling.direct_sparse_discovery import sparse_columns


def commutator(h,word):
    return add(product(h,mono(word)),scale(product(mono(word),h),-1))


def diagonal_commutator(p,energies):
    result={}
    for w,c in p.items():
        value=c*sum((2*flag-1)*energies[j] for flag,j in w)
        if value:result[w]=value
    return result


def polynomial_groups(h,m,enrich=True,creator_channels=False,one_body_steps=0,diagonal_driver=False):
    if type(one_body_steps) is not int or not 0<=one_body_steps<=4:raise ValueError('Krylov depth must be an integer from zero to four')
    if one_body_steps and not enrich:raise ValueError('Krylov enrichment requires commutator generators')
    baseline,signature,_=partition(h,m,'quadratic',True)
    groups=[{'name':g['name'],'polynomials':[canonical(mono(w)) for w in g['words']]} for g in baseline]
    h1={w:c for w,c in h.items() if len(w)==2}
    energies={i:h1.get(((1,i),(0,i)),F(0)) for i in range(m)}
    if enrich:
        for i in range(m):
            whole=commutator(h,((0,i),))
            if any(len(w)>3 for w in whole):raise AssertionError('Two-body commutator exceeded cubic degree')
            cubic={w:c for w,c in whole.items() if len(w)==3}
            if not cubic:continue
            pieces=([ {w:c for w,c in cubic.items() if w[0][1]==creator}
                       for creator in sorted({w[0][1] for w in cubic}) ] if creator_channels else [cubic])
            frontier=list(pieces)
            for _ in range(one_body_steps):
                if diagonal_driver:
                    # [sum_i e_i n_i, w] = sum_(c,i in w) (2c-1)e_i * w.
                    frontier=[diagonal_commutator(p,energies) for p in frontier]
                else:
                    frontier=[{w:c for w,c in add(product(h1,p),scale(product(p,h1),-1)).items() if len(w)==3} for p in frontier]
                frontier=[p for p in frontier if p];pieces.extend(frontier)
            for p in [q for piece in pieces for q in (piece,canonical(adj(piece)))]:
                keys={(sum(2*c-1 for c,_ in w),signature(w)) for w in p}
                if len(keys)!=1:raise ValueError('Inhomogeneous commutator symmetry/charge')
                charge,sig=next(iter(keys))
                target=next(g for g in groups if g['name'].startswith('linear'+('-' if charge<0 else '+'))
                            and signature(next(iter(g['polynomials'][0])))==sig)
                # Power-of-two scaling controls magnitudes without introducing
                # unrelated large prime denominators in the exact exporter.
                largest=max(abs(c) for c in p.values());exponent=float(largest).hex().split('p')[1]
                target['polynomials'].append(scale(p,F(2)**(-int(exponent))))
    return groups,signature


def gram_columns(polys):
    """Upper-triangle Q coefficients; every CAR monomial retained exactly."""
    columns=[];indices=[];word_pair_work=0
    for i,p in enumerate(polys):
        pa=canonical(adj(p))
        for j in range(i,len(polys)):
            q=polys[j];z=product(pa,q);word_pair_work+=len(pa)*len(q)
            if i!=j:z=add(z,canonical(adj(z)))
            if not hermitian(z):raise AssertionError('Real symmetric Gram coefficient is not Hermitian')
            columns.append(z);indices.append(i*len(polys)+j)
    return columns,indices,word_pair_work


def export_groups(groups,grams,rounding=10**9,transforms=None):
    expanded=[];denominator=1;clipped_negative=0.
    if transforms is None:transforms=[None]*len(groups)
    if len(transforms)!=len(groups):raise ValueError("Transform count mismatch")
    for g,q,W in zip(groups,grams,transforms):
        ev,u=np.linalg.eigh((q+q.T)/2);clipped_negative+=float(-ev[ev<0].sum())
        effective=np.sqrt(np.maximum(ev,0))[:,None]*u.T
        if W is not None:
            W=np.asarray(W,dtype=float)
            if W.shape!=q.shape or not np.all(np.isfinite(W)):raise ValueError('Invalid export transform')
            effective=effective@W
        if not np.all(np.isfinite(effective)) or np.max(np.abs(effective),initial=0)*rounding>=2**62:
            raise ValueError('Nonfinite or excessive export factor')
        z=np.rint(effective*rounding).astype(np.int64)
        rows=[]
        for row in z:
            p=add(*(scale(poly,F(int(coef),rounding)) for coef,poly in zip(row,g['polynomials']) if coef))
            if p:rows.append(p)
            for c in p.values():denominator=lcm(denominator,c.denominator)
        expanded.append((g['name'],rows))
    blocks=[]
    for name,rows in expanded:
        if not rows:continue
        words=sorted({w for p in rows for w in p},key=lambda w:(len(w),w))
        if len({sum(2*c-1 for c,_ in w) for w in words})!=1:raise AssertionError('Export charge mixed')
        factors=[[int(p.get(w,F(0))*denominator) for w in words] for p in rows]
        blocks.append({'name':name,'words':words,'factor':factors})
    return blocks,denominator,clipped_negative


def run(h,m,n,out,solver_seconds=60,enrich=True,creator_channels=False,one_body_steps=0,diagonal_driver=False,solver='CLARABEL',map_backend='exact',row_condition=False,basis_condition='none'):
    if solver not in ('CLARABEL','SCS'):raise ValueError('Unsupported solver')
    if map_backend not in ('exact','contraction'):raise ValueError('Unsupported map backend')
    if basis_condition not in ('none','whiten') or (basis_condition=='whiten' and map_backend!='contraction'):raise ValueError('Whitening requires the contraction backend')
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    if not hermitian(h) or any(len(w)>4 for w in h):raise ValueError('Hermitian two-body H required')
    groups,signature=polynomial_groups(h,m,enrich,creator_channels,one_body_steps,diagonal_driver);zero=signature(())
    transforms=[None]*len(groups);conditioning=[]
    if basis_condition=='whiten':
        from research.certificate_scaling.polynomial_gram_contraction import whitening_transform
        for i,g in enumerate(groups):
            W,stats=whitening_transform(g['polynomials']);transforms[i]=W
            conditioning.append({'name':g['name'],**stats})
        (out/'basis_transforms.json').write_text(json.dumps([W.tolist() for W in transforms])+'\n')
    basis=[p for p in multiplier_basis(m,max_body=2) if all(signature(w)==zero for w in p)]
    free_polys=[mono(())]+[product(number_shift(m,n),p) for p in basis]
    columns=[];mapspec=[];work=0;map_stats={}
    if map_backend=='exact':
        for g in groups:
            cs,ix,cost=gram_columns(g['polynomials']);columns.append(cs);mapspec.append(ix);work+=cost
        mapped_words={w for cs in columns for p in cs for w in p}
    else:
        from research.certificate_scaling.polynomial_gram_contraction import prepare,contract
        prepared,mapped_words,map_stats=prepare(groups);work=map_stats['monomial_word_pairs']
    allwords=set(h)|{w for p in free_polys for w in p}|mapped_words
    if any(len(w)>6 or sum(2*c-1 for c,_ in w)!=0 for w in allwords):raise AssertionError('Unexpected square degree or charge')
    rows=sorted(allwords,key=lambda w:(len(w),w));lookup={w:i for i,w in enumerate(rows)}
    free=sparse_columns(free_polys,lookup);rhs=np.array([float(h.get(w,0)) for w in rows])
    x=cp.Variable(len(free_polys));rem=cp.Variable(len(rows));expr=free@x;variables=[];nnz=0
    squared=np.asarray(free.power(2).sum(axis=1)).ravel()
    for gid,g in enumerate(groups):
        k=len(g['polynomials']);q=cp.Variable((k,k),PSD=True);variables.append(q)
        if map_backend=='exact':matrix=sparse_columns(columns[gid],lookup);ix=mapspec[gid]
        else:
            matrix,ix,stats=contract(prepared[gid],g['polynomials'],lookup,transform=transforms[gid])
            map_stats.setdefault('blocks',[]).append(stats)
        nnz+=matrix.nnz;expr+=matrix@cp.reshape(q,(k*k,),order='C')[ix]
        if row_condition:squared+=np.asarray(matrix.power(2).sum(axis=1)).ravel()
    row_scale=1/np.maximum(1.,np.sqrt(squared)) if row_condition else np.ones(len(rows))
    problem=cp.Problem(cp.Minimize(cp.norm1(rem)-x[0]),[cp.multiply(row_scale,expr+rem-rhs)==0,rem[lookup[()]]==0])
    built=time.monotonic()
    (out/'pre_solve.json').write_text(json.dumps({'stage':'solver_starting','construction_seconds':built-start,
        'modes':m,'particles':n,'one_body_steps':one_body_steps,'diagonal_driver':diagonal_driver,
        'polynomial_Gram_entries':sum(len(g['polynomials'])**2 for g in groups),
        'symbolic_word_pair_products':work,'coefficient_rows':len(rows),'gram_map_nonzeros':nnz,
        'basis_condition':basis_condition,'basis_conditioning':conditioning,'solver':solver,'map_backend':map_backend,'map_stats':map_stats,'row_condition':row_condition,'solver_seconds_budget':solver_seconds,'scope':'Construction receipt only; no completed solve or accepted certificate claimed.'},indent=2)+'\n')
    options=({'tol_gap_abs':1e-8,'tol_feas':1e-8,'tol_gap_rel':1e-8,'max_iter':1000,'time_limit':float(solver_seconds)}
             if solver=='CLARABEL' else {'eps':1e-8,'max_iters':150000,'time_limit_secs':float(solver_seconds)})
    problem.solve(solver=solver,**options)
    solved=time.monotonic()
    if x.value is None or any(q.value is None for q in variables):raise RuntimeError('No bounded Gram solution')
    if problem.constraints[0].dual_value is not None:
        (out/'dual_proposal.json').write_text(json.dumps({'modes':m,'particles':n,
            'hamiltonian':encode(h),'enrich':enrich,'creator_channels':creator_channels,'one_body_steps':one_body_steps,'diagonal_driver':diagonal_driver,
            'rows':[[list(letter) for letter in w] for w in rows],
            'dual':[float(v) for v in problem.constraints[0].dual_value*row_scale],
            'scope':'Floating proposal only; not an accepted dual witness.'},separators=(',',':'))+'\n')
    blocks,den,negative=export_groups(groups,[q.value for q in variables],transforms=transforms);X=add(*(scale(p,F(round(float(v)*10**12),10**12)) for p,v in zip(basis,x.value[1:])))
    cert={'modes':m,'particles':n,'hamiltonian':encode(h),'b':str(F(round(float(x.value[0])*10**12),10**12)),'number_multiplier':encode(X),'denominator':den,'blocks':blocks}
    exact=verify(cert);raw=json.dumps(cert,separators=(',',':'))+'\n';(out/'certificate.json').write_text(raw)
    receipt={'method':'H_only_cubic_commutator_plus_quadratic' if enrich else 'quadratic_baseline_body2_ideal',
             'modes':m,'particles':n,'polynomial_Gram_dimensions':[len(g['polynomials']) for g in groups],
             'polynomial_Gram_entries':sum(len(g['polynomials'])**2 for g in groups),
             'polynomial_term_counts':[[len(p) for p in g['polynomials']] for g in groups],
             'symbolic_word_pair_products':work,'coefficient_rows':len(rows),'coefficient_max_degree':max(map(len,rows)),
             'gram_map_nonzeros':nnz,'map_backend':map_backend,'map_stats':map_stats,'ideal_variables':len(basis),'construction_seconds':built-start,
             'solve_seconds':solved-built,'export_seconds':time.monotonic()-solved,'wall_seconds':time.monotonic()-start,
             'status':problem.status,'solver':solver,'solver_options':options,'row_condition':row_condition,'numeric_lower':-float(problem.value),'raw_b':float(x.value[0]),'negative_eigenvalue_mass':negative,
             'exact':exact,'certificate_bytes':len(raw.encode()),'factor_denominator_bits':den.bit_length(),
             'basis_condition':basis_condition,'basis_conditioning':conditioning,'creator_channels':creator_channels,'one_body_steps':one_body_steps,'diagonal_driver':diagonal_driver,'source_factors_used':False,'source_upper_used':False,'scope':'Fixed H-only polynomial span; no omitted-cone optimality or generic accuracy guarantee.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k not in ('polynomial_term_counts','polynomial_Gram_dimensions')}),flush=True);return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--solver-seconds',type=float,default=60);p.add_argument('--baseline',action='store_true');p.add_argument('--creator-channels',action='store_true');p.add_argument('--one-body-steps',type=int,default=0);p.add_argument('--diagonal-driver',action='store_true');p.add_argument('--solver',choices=['CLARABEL','SCS'],default='CLARABEL');p.add_argument('--map-backend',choices=['exact','contraction'],default='exact');p.add_argument('--row-condition',action='store_true');p.add_argument('--basis-condition',choices=['none','whiten'],default='none');a=p.parse_args();f=json.loads(a.fixture.read_text());run(decode(f['hamiltonian'],f['modes'],4),f['modes'],f['particles'],a.out,a.solver_seconds,not a.baseline,a.creator_channels,a.one_body_steps,a.diagonal_driver,a.solver,a.map_backend,a.row_condition,a.basis_condition)
