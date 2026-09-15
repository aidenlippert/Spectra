"""Energy-scored quartic column generation over a symmetry-reduced cubic SDP.

No full quartic PSD optimization supplies the square directions. Candidate
directions come from the current restricted problem's moment dual; small
lookahead optimizations score their actual effect on the energy objective.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import time

import cvxpy as cp
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

from experiments.marginal_coefficient import gram_map, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, upper_bound
from experiments.marginal_dual_exact import component_charge, full_word_blocks
from experiments.marginal_reynolds import invariant_rows, reynolds_matrix
from experiments.marginal_stabilizer import split_blocks
from experiments.marginal_orbit_quartic import symmetry_group
from experiments.marginal_orbit_certificate import normalize_scalar_residual, average_canonical_polynomial
from experiments.marginal_symbolic import add, encode, multiplier_basis, number_shift, product, scale


def representatives(blocks,modes):
    seen=set();result=[]
    for i,block in enumerate(blocks):
        key=tuple(sorted(component_charge(block['words'][0],modes)))
        if key not in seen:seen.add(key);result.append(i)
    return result


def prepare(modes,pricing='full'):
    if pricing not in ('full','stabilizer'):raise ValueError('Invalid pricing method')
    h=hopping_polynomial(modes,F(1,5));rows=invariant_rows(modes,4)
    projection,chosen=reynolds_matrix(rows,modes);reduced=projection[chosen,:]
    lookup={w:i for i,w in enumerate(rows)}
    basis=[p for p in multiplier_basis(modes,3) if all(not any(component_charge(w,modes)) for w in p)]
    ideal=[product(number_shift(modes,modes//2),p) for p in basis]
    a=reduced@np.array([[float(col.get(w,0)) for col in ideal] for w in rows])
    _,r,piv=la.qr(a,mode='economic',pivoting=True);rank=int(np.sum(abs(np.diag(r))>1e-10));piv=piv[:rank]
    cubic=full_word_blocks(modes,3);indices=set(representatives(cubic,modes));split,_=split_blocks(modes,3)
    base=[b for b in split if int(b['name'].split(':')[0]) in indices];base_maps=[]
    for block in base:
        tr=sp.csr_matrix(block['basis_transform'])
        base_maps.append(reduced@gram_map(block['words'],lookup)@sp.kron(tr,tr,format='csr'))
    quartic=full_word_blocks(modes,4)
    candidates=[quartic[i] for i in representatives(quartic,modes)]
    maps=[reduced@gram_map(b['words'],lookup) for b in candidates]
    pricing_blocks=[]
    if pricing=='stabilizer':
        candidate_indices=representatives(quartic,modes)
        positions={index:i for i,index in enumerate(candidate_indices)}
        pieces,_=split_blocks(modes,4,indices=candidate_indices)
        for piece in pieces:
            index=positions[int(piece['name'].split(':')[0])]
            tr=sp.csr_matrix(piece['basis_transform'])
            pricing_blocks.append({'candidate':index,'transform':piece['basis_transform'],
                                   'map':maps[index]@sp.kron(tr,tr,format='csr')})
    pricing_groups={}
    for block in pricing_blocks:
        pricing_groups.setdefault(block['transform'].shape[1],[]).append(block)
    pricing_groups=[{'dimension':d,'blocks':blocks,'map':sp.hstack([b['map'] for b in blocks],format='csr')}
                    for d,blocks in sorted(pricing_groups.items())]
    weights=[]
    for i in chosen:
        weights.append(sum(1 if tuple(j for c,j in rows[k] if c)==tuple(j for c,j in rows[k] if not c) else 2 for k in projection.getrow(i).indices))
    return {'modes':modes,'h':h,'rows':rows,'reduced':reduced,'base':base,'base_maps':base_maps,
            'candidates':candidates,'candidate_maps':maps,'a':a[:,piv],
            'basis':[basis[int(j)] for j in piv],'weights':np.array(weights),
            'rhs':reduced@np.array([float(h.get(w,0)) for w in rows]),
            'unit':reduced@np.array([float(w==()) for w in rows]),
            'pricing':pricing,'pricing_blocks':pricing_blocks,'pricing_groups':pricing_groups}


def restricted_problem(data,cuts,pool_size):
    n=len(data['rhs']);b=cp.Variable();x=cp.Variable(len(data['basis']));residual=cp.Variable(n)
    dimensions=[block['basis_transform'].shape[1] for block in data['base']]
    grams=[cp.Variable((d,d),PSD=True) for d in dimensions]
    expression=b*data['unit']+data['a']@x+residual
    expression+=sp.hstack(data['base_maps'],format='csr')@cp.hstack([cp.reshape(q,(d*d,),order='C') for q,d in zip(grams,dimensions)])
    cut_weights=cp.Variable(len(cuts),nonneg=True) if cuts else None
    extra=[]
    if cuts:
        expression+=np.column_stack([c['column'] for c in cuts])@cut_weights;extra.append(cut_weights<=100)
    trial=cp.Parameter((n,pool_size),value=np.zeros((n,pool_size)))
    trial_weights=cp.Variable(pool_size,nonneg=True);expression+=trial@trial_weights;extra.append(trial_weights<=100)
    equality=expression==data['rhs']
    problem=cp.Problem(cp.Maximize(b-cp.norm1(cp.multiply(data['weights'],residual))),[equality,*extra])
    def solve(columns=None,eps=1e-9):
        trial.value=np.zeros((n,pool_size))
        if columns:
            matrix=np.zeros((n,pool_size));matrix[:,:len(columns)]=np.column_stack(columns);trial.value=matrix
        problem.solve(solver='CLARABEL',tol_gap_abs=eps,tol_gap_rel=eps,tol_feas=eps,max_iter=150)
        if b.value is None or equality.dual_value is None or any(q.value is None for q in grams):raise RuntimeError(problem.status)
        return {'objective':float(problem.value),'b':float(b.value),'x':np.array(x.value),
                'grams':[np.array(q.value) for q in grams],'cut_weights':np.array(cut_weights.value) if cuts else np.array([]),
                'dual':np.array(equality.dual_value),'status':problem.status}
    return solve


def price_candidates(data,dual,pool_size):
    proposals=[]
    if data.get('pricing')=='stabilizer':
        best={}
        for group in data['pricing_groups']:
            d=group['dimension'];matrices=np.asarray(group['map'].T@dual).reshape(-1,d,d)
            values,vectors=np.linalg.eigh((matrices+matrices.transpose(0,2,1))/2)
            for j,block in enumerate(group['blocks']):
                i=block['candidate'];value=float(values[j,0])
                if value<-1e-7 and (i not in best or value<best[i]['violation']):
                    best[i]={'candidate':i,'violation':value,'small_vector':vectors[j,:,0],'block':block}
        for p in sorted(best.values(),key=lambda p:p['violation'])[:pool_size]:
            v=p.pop('small_vector');block=p.pop('block')
            p['vector']=block['transform']@v
            p['column']=np.asarray(block['map']@np.outer(v,v).ravel()).ravel()
            proposals.append(p)
        return proposals,len(best)
    for i,(block,g) in enumerate(zip(data['candidates'],data['candidate_maps'])):
        n=len(block['words']);matrix=np.asarray(g.T@dual).reshape(n,n)
        values,vectors=np.linalg.eigh((matrix+matrix.T)/2)
        if values[0]<-1e-7:
            v=vectors[:,0];column=np.asarray(g@np.outer(v,v).ravel()).ravel()
            proposals.append({'candidate':i,'violation':float(values[0]),'vector':v,'column':column})
    proposals.sort(key=lambda p:p['violation'])
    return proposals[:pool_size],len(proposals)


def exact_certificate(data,cuts,solution):
    modes=data['modes'];squares=[]
    def add_square(words,vector,weight):
        w=F(round(max(0.,float(weight))*10**12),10**12)
        p={word:F(round(float(c)*10**12),10**12) for word,c in zip(words,vector) if round(float(c)*10**12)}
        if p and w:squares.append({'polynomial':encode(p),'weight':str(w)})
    for block,q in zip(data['base'],solution['grams']):
        values,vectors=np.linalg.eigh((q+q.T)/2)
        for j,value in enumerate(values):
            if value>0:add_square(block['words'],block['basis_transform']@vectors[:,j]*np.sqrt(value),1)
    for cut,weight in zip(cuts,solution['cut_weights']):
        add_square(data['candidates'][cut['candidate']]['words'],cut['vector'],weight)
    raw_x=add(*(scale(p,F(round(float(c)*10**12),10**12)) for p,c in zip(data['basis'],solution['x'])))
    group=symmetry_group(modes);x=average_canonical_polynomial(raw_x,group)
    certificate={'modes':modes,'particles':modes//2,'operator_degree':4,'hamiltonian':encode(data['h']),
                 'number_multiplier':encode(x),'b':str(F(round(solution['b']*10**12),10**12)),
                 'permutations':[list(g) for g in group],'orbit_squares':squares,
                 'variational_upper':symmetric_upper(modes,F(1,5))}
    certificate,receipt=normalize_scalar_residual(certificate)
    upper=upper_bound(modes,F(1,5),certificate['variational_upper']['amplitudes'])['upper']
    receipt.update({'upper':upper,'width':str(F(upper)-F(receipt['lower'])),'numerical_objective':solution['objective']})
    receipt['width_float']=float(F(receipt['width']))
    if F(receipt['width'])<0:raise ValueError('Inconsistent interval')
    return certificate,receipt


def run(modes=10,iterations=20,pool_size=8,resume=None,strategy='energy',max_cuts=None,pricing='full'):
    if strategy not in ('energy','negative','batch'):raise ValueError('Invalid selection strategy')
    if type(iterations) is not int or iterations<0 or type(pool_size) is not int or pool_size<1:raise ValueError('Invalid search budget')
    if max_cuts is not None and (type(max_cuts) is not int or max_cuts<0):raise ValueError('Invalid cut budget')
    started=time.monotonic();data=prepare(modes,pricing);cuts=[];history=[]
    start_iteration=0
    if resume is not None:
        state=json.loads(Path(resume).read_text())
        if state['modes']!=modes or state['t']!='1/5':raise ValueError('Checkpoint model mismatch')
        start_iteration=state['iteration']
        if not 0<=start_iteration<=iterations:raise ValueError('Checkpoint iteration outside run budget')
        for saved in state['cuts']:
            index=saved['candidate'];v=np.array(saved['vector'],dtype=float)
            if v.shape!=(len(data['candidates'][index]['words']),) or not np.all(np.isfinite(v)):raise ValueError('Invalid checkpoint direction')
            cuts.append({'candidate':index,'vector':v,'column':np.asarray(data['candidate_maps'][index]@np.outer(v,v).ravel()).ravel()})
    suffix='' if strategy=='energy' and max_cuts is None else f'_{strategy}_cap{max_cuts}'
    if pricing!='full':suffix+=f'_{pricing}'
    out=Path(__file__).resolve().parents[1]/'results/marginal_energy_selector'/f'm{modes}_steps{iterations}{suffix}'
    out.mkdir(parents=True,exist_ok=True)
    structure={'modes':modes,'coefficient_rows':len(data['rhs']),'cubic_psd_blocks':len(data['base']),
               'largest_cubic_psd_block':max(b['basis_transform'].shape[1] for b in data['base']),
               'cubic_psd_scalars':sum(b['basis_transform'].shape[1]*(b['basis_transform'].shape[1]+1)//2 for b in data['base']),
               'quartic_candidate_dictionaries':len(data['candidates']),'build_seconds':time.monotonic()-started,
               'strategy':strategy,'maximum_cut_weight':100,'max_cuts':max_cuts,
               'pricing':pricing,'pricing_psd_blocks':len(data['pricing_blocks']) if pricing=='stabilizer' else len(data['candidates']),
               'largest_pricing_matrix':max(b['transform'].shape[1] for b in data['pricing_blocks']) if pricing=='stabilizer' else max(len(b['words']) for b in data['candidates'])}
    (out/'structure.json').write_text(json.dumps(structure,indent=2)+'\n');print(json.dumps(structure),flush=True)
    def checkpoint(iteration):
        state={'modes':modes,'t':'1/5','iteration':iteration,
               'cuts':[{'candidate':p['candidate'],'vector':p['vector'].tolist()} for p in cuts]}
        (out/'cut_state.json').write_text(json.dumps(state)+'\n')
    for iteration in range(start_iteration,iterations+1):
        start=time.monotonic();solve=restricted_problem(data,cuts,pool_size);base=solve()
        pool,violated=price_candidates(data,base['dual'],pool_size)
        step={'iteration':iteration,'cuts':len(cuts),'objective':base['objective'],'status':base['status'],
              'violated_dictionaries':violated,'worst_eigenvalue':pool[0]['violation'] if pool else 0}
        if iteration==start_iteration:
            cert,receipt=exact_certificate(data,cuts,base)
            (out/'initial_certificate.json').write_text(json.dumps(cert)+'\n');(out/'initial_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if iteration==iterations or not pool or (max_cuts is not None and len(cuts)>=max_cuts):
            step['stop']='no_detected_violation' if not pool else ('cut_budget' if max_cuts is not None and len(cuts)>=max_cuts else 'iteration_limit')
            history.append(step);break
        scores=[]
        if strategy=='negative':selected=[pool[0]];step['selection']='most_negative'
        elif strategy=='batch':selected=pool;step['selection']='most_negative_batch'
        else:
            for candidate in pool:
                result=solve([candidate['column']]);scores.append(result['objective'])
            best=int(np.argmax(scores));gain=scores[best]-base['objective']
            if gain>1e-7:
                selected=[pool[best]];step['selection']='best_single'
            else:
                batch=solve([p['column'] for p in pool]);step['batch_objective']=batch['objective']
                selected=pool;step['selection']='batch_after_flat_singles'
        if max_cuts is not None:selected=selected[:max_cuts-len(cuts)]
        step.update({'single_gains':[score-base['objective'] for score in scores],
                     'selected_candidates':[p['candidate'] for p in selected],'seconds':time.monotonic()-start})
        cuts.extend(selected);history.append(step)
        checkpoint(iteration+1)
        (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        print(json.dumps({k:step[k] for k in ('iteration','cuts','objective','selection','seconds')}),flush=True)
    certificate,receipt=exact_certificate(data,cuts,base)
    checkpoint(iteration)
    receipt.update({'cuts':len(cuts),'iterations':len(history)-1,'last_iteration':iteration,'resumed_from':str(resume) if resume else None,'total_seconds':time.monotonic()-started,
                    'strategy':strategy,'maximum_cut_weight':100,'max_cuts':max_cuts,'pricing':pricing})
    (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    (out/'certificate.json').write_text(json.dumps(certificate)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
    return certificate,receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--modes',type=int,default=10);parser.add_argument('--iterations',type=int,default=20);parser.add_argument('--pool',type=int,default=8);parser.add_argument('--resume',type=Path);parser.add_argument('--strategy',choices=['energy','negative','batch'],default='energy');parser.add_argument('--max-cuts',type=int)
    parser.add_argument('--pricing',choices=['full','stabilizer'],default='full')
    args=parser.parse_args();run(args.modes,args.iterations,args.pool,args.resume,args.strategy,args.max_cuts,args.pricing)
