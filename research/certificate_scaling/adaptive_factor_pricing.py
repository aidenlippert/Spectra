"""Dual-guided sparse CAR square discovery with explicit full pricing costs.

No solved Gram or wavefunction enters discovery. Pricing constructs complete
moment maps and numerical moment matrices: this is NOT an omitted-family
complexity theorem. Every accepted energy floor is independently rationalized.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from experiments.marginal_coefficient import dictionaries, gram_map, dagger
from experiments.marginal_symbolic import add, canonical, decode, encode, hermitian, mono, multiplier_basis, number_shift, product, scale, verify, word_product
from research.certificate_scaling.direct_sparse_discovery import candidates, sparse_columns


def pricing_vectors(matrix, width, count):
    """Best pair seeds; greedy larger support. Width=0 uses full eigenvectors."""
    n=len(matrix); out=[]
    if width==0:
        vals,vec=np.linalg.eigh(matrix)
        for i in np.argsort(vals)[:count]:
            if vals[i]<-1e-8:out.append((float(vals[i]),np.arange(n),vec[:,i]))
        return out
    width=min(width,n)
    if n==1:
        return [(float(matrix[0,0]),np.array([0]),np.array([1.]))] if matrix[0,0]<-1e-8 else []
    ii,jj=np.triu_indices(n,1); di=np.diag(matrix)
    value=(di[ii]+di[jj]-np.sqrt((di[ii]-di[jj])**2+4*matrix[ii,jj]**2))/2
    seeds=np.argsort(value)[:count]
    for seed in seeds:
        support=[int(ii[seed]),int(jj[seed])]
        for step in range(max(1,width-1)):
            vals,vec=np.linalg.eigh(matrix[np.ix_(support,support)])
            if len(support)>=width:break
            mask=np.ones(n,dtype=bool);mask[support]=False
            available=np.flatnonzero(mask)
            coupling=matrix[np.ix_(available,support)]@vec[:,0]
            score=coupling**2/np.maximum(di[available]-vals[0],1e-8)
            support.append(int(available[np.argmax(score)]))
        vals,vec=np.linalg.eigh(matrix[np.ix_(support,support)])
        if vals[0]<-1e-8:out.append((float(vals[0]),np.array(support),vec[:,0]))
    return out


def run(h,modes,particles,outdir,family='quadratic',width=2,iterations=20,
        batch=32,atom_cap=4096,ideal_body=1,seconds=120,denominator=10**7,
        half_rows=False,prune=False):
    started=time.monotonic();h=canonical(h)
    if not hermitian(h) or any(len(w)>4 or sum(2*c-1 for c,_ in w) for w in h):
        raise ValueError('Hermitian, number-conserving two-body H required')
    if width not in (0,2,4,8,16,32) or not 1<=iterations<=200 or atom_cap<512:
        raise ValueError('Invalid pricing budget')
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    blocks=dictionaries(modes,family)
    max_body=2 if family=='quadratic' and ideal_body<=1 else 3
    rows=[tuple((1,p) for p in left)+tuple((0,p) for p in right)
          for k in range(max_body+1) for left in combinations(range(modes),k)
          for right in combinations(range(modes),k) if not half_rows or left<=right]
    residual_weights=np.array([1 if not half_rows or tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in rows])
    lookup={w:i for i,w in enumerate(rows)}
    maps=[gram_map(b['words'],lookup) for b in blocks]
    basis=multiplier_basis(modes,max_body=ideal_body)
    shift=number_shift(modes,particles)
    ideals=[product(shift,q) for q in basis]
    free=sparse_columns([{w:c for w,c in p.items() if w in lookup} for p in [mono(())]+ideals],lookup)
    rhs=np.array([float(h.get(w,0)) for w in rows]); nfree=free.shape[1]
    atoms=[];columns=[];seen=set();zero_cross_added=0;pruned_count=0;cumulative_atoms=0
    def atom_key(gid,support,integers):
        div=math.gcd(*map(abs,integers));sgn=1 if integers[0]>0 else -1
        return (gid,tuple(support),tuple(v//div*sgn for v in integers))
    def append_atom(gid,support,integers,dd=10000):
        nonlocal cumulative_atoms
        nonzero=[(int(i),int(v)) for i,v in zip(support,integers) if int(v)]
        if not nonzero:return False
        nonzero.sort();support=np.array([x[0] for x in nonzero]);integers=np.array([x[1] for x in nonzero])
        key=atom_key(gid,support.tolist(),integers.tolist())
        if key in seen:return False
        n=len(blocks[gid]['words']);coeff=integers.astype(float)/dd
        indices=(support[:,None]*n+support[None,:]).ravel()
        col=maps[gid][:,indices]@sparse.csc_matrix(np.outer(coeff,coeff).ravel()[:,None])
        atoms.append((gid,support.tolist(),integers.tolist(),dd));columns.append(col);seen.add(key);cumulative_atoms+=1
        return True
    seed,seed_counts=candidates(h,modes,512)
    for ws,signs,_ in seed:
        for gid,block in enumerate(blocks):
            if all(w in block['words'] for w in ws):
                append_atom(gid,[block['words'].index(w) for w in ws],[s*10000 for s in signs]);break
    map_seconds=time.monotonic()-started;history=[];best=None;bestrec=None
    total_solve=0.;total_price=0.;total_export=0.;stop='iteration_budget'
    for iteration in range(iterations+1):
        if time.monotonic()-started>seconds and best is not None:
            stop='wall_budget';break
        before=time.monotonic()
        A=sparse.hstack([free]+columns,format='csc')
        eq=sparse.hstack([A,sparse.eye(len(rows)),-sparse.eye(len(rows))],format='csc')
        objective=np.r_[-1.,np.zeros(nfree-1+len(atoms)),residual_weights,residual_weights]
        result=linprog(objective,A_eq=eq,b_eq=rhs,
            bounds=[(None,None)]*nfree+[(0,None)]*(len(atoms)+2*len(rows)),
            method='highs',options={'time_limit':max(1.,min(60.,seconds-(time.monotonic()-started)))})
        total_solve+=time.monotonic()-before
        if not result.success:
            if best is None:raise RuntimeError('Initial LP failed: '+result.message)
            stop='LP_limit_or_failure';break
        # Min-LP marginal has constant entry -1; physical y is its negative.
        y=-np.asarray(result.eqlin.marginals)
        row={'iteration':iteration,'atoms':len(atoms),'numeric_lower':float(-result.fun),
             'LP_iterations':int(result.nit),'LP_nonzeros':int(eq.nnz)}
        before=time.monotonic()
        fs=[]
        for weight,(gid,support,ints,dd) in zip(result.x[nfree:A.shape[1]],atoms):
            factor=[int(round(math.sqrt(max(0,weight))*v/dd*denominator)) for v in ints]
            if any(factor):fs.append({'name':f'priced-{len(fs)}','words':[blocks[gid]['words'][i] for i in support],'factor':[factor]})
        X=add(*(scale(q,F(int(round(v*10**9)),10**9)) for q,v in zip(basis,result.x[1:nfree])))
        cert={'modes':modes,'particles':particles,'hamiltonian':encode(h),
              'b':str(F(int(round(result.x[0]*10**9)),10**9)),'number_multiplier':encode(X),
              'denominator':denominator,'blocks':fs}
        rec=verify(cert);row['exact_lower']=rec['lower_float']
        if best is None or F(rec['lower'])>F(bestrec['lower']):best,bestrec=cert,rec
        total_export+=time.monotonic()-before
        history.append(row)
        (outdir/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        (outdir/'certificate.json').write_text(json.dumps(best,separators=(',',':'))+'\n')
        if iteration==iterations:break
        if prune and iteration%5==0:
            keep=[]
            for j,(atom,weight) in enumerate(zip(atoms,result.x[nfree:A.shape[1]])):
                if weight>1e-9:keep.append(j)
                else:
                    seen.remove(atom_key(atom[0],atom[1],atom[2]));pruned_count+=1
            atoms=[atoms[j] for j in keep];columns=[columns[j] for j in keep]
        if len(atoms)>=atom_cap:stop='atom_budget';break
        before=time.monotonic();proposals=[]
        for gid,gmap in enumerate(maps):
            n=len(blocks[gid]['words']);mat=np.asarray(gmap.T@y).reshape(n,n);mat=(mat+mat.T)/2
            for val,support,vec in pricing_vectors(mat,width,batch):
                proposals.append((val,gid,support,np.rint(vec*10000).astype(int)))
        added=0;zero_cross=0
        for val,gid,support,ints in sorted(proposals,key=lambda p:p[0]):
            if len(atoms)>=atom_cap or added>=batch:break
            if append_atom(gid,support,ints):
                added+=1
                if len(support)==2:
                    u,v=[blocks[gid]['words'][i] for i in support]
                    cross=add(dict(word_product(dagger(u),v)),dict(word_product(dagger(v),u)))
                    if sum(c*h.get(w,F(0)) for w,c in cross.items())==0:zero_cross+=1
        zero_cross_added+=zero_cross;total_price+=time.monotonic()-before
        row.update({'pricing_proposals':len(proposals),'added_atoms':added,'added_zero_H_cross_pairs':zero_cross,
                    'minimum_proposed_dual_value':min((p[0] for p in proposals),default=None)})
        if not added:stop='no_new_heuristic_direction';break
    bestrec.update({'method':'adaptive_dual_factor_pricing','family':family,'width':width,
        'ideal_body':ideal_body,'atom_cap':atom_cap,'atoms_retained':len(atoms),'stop':stop,
        'Hermitian_half_rows':half_rows,'prune_inactive_atoms':prune,
        'cumulative_atoms_added':cumulative_atoms,'cumulative_atoms_pruned':pruned_count,
        'rounds_solved':len(history),'seed_counts':seed_counts,'coefficient_rows':len(rows),
        'complete_pricing_Gram_dimensions':[len(b['words']) for b in blocks],
        'complete_pricing_matrix_entries':sum(len(b['words'])**2 for b in blocks),
        'complete_pricing_map_nonzeros':sum(m.nnz for m in maps),
        'added_zero_H_cross_pairs':zero_cross_added,'map_build_seconds':map_seconds,
        'total_LP_seconds':total_solve,'total_pricing_seconds':total_price,
        'total_exact_export_seconds':total_export,'wall_seconds':time.monotonic()-started,
        'source_factors_used':False,'source_upper_used':False,'many_body_space_enumerated':False,
        'omitted_family_optimality_proved':False,
        'certificate_bytes':len(json.dumps(best,separators=(',',':')).encode())})
    (outdir/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    (outdir/'receipt.json').write_text(json.dumps(bestrec,indent=2)+'\n')
    return best,bestrec


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True);p.add_argument('--outputdir',type=Path,required=True)
    p.add_argument('--family',choices=['quadratic','local','mixed'],default='quadratic')
    p.add_argument('--width',type=int,default=2);p.add_argument('--iterations',type=int,default=20)
    p.add_argument('--batch',type=int,default=32);p.add_argument('--atom-cap',type=int,default=4096)
    p.add_argument('--ideal-body',type=int,default=1);p.add_argument('--seconds',type=float,default=120)
    p.add_argument('--half-rows',action='store_true');p.add_argument('--prune',action='store_true')
    a=p.parse_args();f=json.loads(a.fixture.read_text())
    _,r=run(decode(f['hamiltonian'],f['modes'],4),f['modes'],f['particles'],a.outputdir,
        a.family,a.width,a.iterations,a.batch,a.atom_cap,a.ideal_body,a.seconds,
        half_rows=a.half_rows,prune=a.prune)
    print(json.dumps(r))
