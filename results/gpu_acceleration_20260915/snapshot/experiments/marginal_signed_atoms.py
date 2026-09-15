"""Exact certificates for signed and rational-amplitude sparse positive atoms.

Unlike the balanced-clique certificate, an atom may have either sign on every
coordinate. The versioned rational format also permits unequal amplitudes on
two to four coordinates. The replay path subtracts the atoms
from ``W H W`` using exact rationals and proves the remainder by diagonal
dominance in the metric ``W^2``.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time

from experiments.marginal_clique_gap import metric_matrix
from experiments.marginal_implicit_certificate import rational_text


def _signs(indices, signs, n):
    if (type(indices) is not list or type(signs) is not list
            or len(indices) not in (3, 4) or len(signs) != len(indices)
            or indices != sorted(set(indices))
            or any(type(i) is not int or not 0 <= i < n for i in indices)
            or any(type(s) is not int or s not in (-1, 1) for s in signs)
            or signs[0] != 1):
        raise ValueError('Sorted distinct support-three/four atoms with +/-1 signs anchored at +1 required')
    return tuple(indices), tuple(signs)


def atom_vector(atom,n,allow_amplitudes=False):
    if 'signs' in atom:
        if 'amplitudes' in atom or 'amplitude_scale' in atom: raise ValueError('One vector encoding required')
        return _signs(atom.get('indices'),atom['signs'],n)
    ind,values,scale=atom.get('indices'),atom.get('amplitudes'),atom.get('amplitude_scale')
    if (not allow_amplitudes or type(ind) is not list or not 2<=len(ind)<=4
            or any(type(i) is not int or not 0<=i<n for i in ind) or ind!=sorted(set(ind))
            or type(values) is not list or len(values)!=len(ind)
            or any(type(v) is not int or not 0<abs(v)<=10**9 for v in values)
            or values[0]<=0 or type(scale) is not int or not 1<=scale<=10**9):
        raise ValueError('Bounded rational two/four-support vector with positive first amplitude required')
    return tuple(ind),tuple(F(v,scale) for v in values)


def encode_vector(key):
    """Canonical direction encoding, shared by proofs and search dictionaries."""
    import math
    group,vector=key
    atom={'indices':list(group)}
    if len(group)>=3 and all(x in (-1,1) for x in vector):
        atom['signs']=[int(x) for x in vector]
    else:
        denominator=math.lcm(*(F(x).denominator for x in vector))
        atom.update(amplitude_scale=denominator,amplitudes=[int(x*denominator) for x in vector])
    return atom


def residual(a, weights, atoms, denominator,allow_amplitudes=False):
    """Return ``R=WHW-sum(atom)`` and the exact metric diagonal."""
    if type(denominator) is not int or not 1 <= denominator <= 10**18:
        raise ValueError('Bounded positive integer atom scale required')
    h, metric = metric_matrix(a, weights)
    n = len(h)
    if type(atoms) is not list or len(atoms) > 100000:
        raise ValueError('Bounded explicit signed atom list required')
    r = [row[:] for row in h]
    seen = set()
    for atom in atoms:
        if type(atom) is not dict:
            raise ValueError('Signed atom object required')
        ind, sig = atom_vector(atom,n,allow_amplitudes)
        w = atom.get('weight')
        if type(w) is not int or w <= 0:
            raise ValueError('Signed atom weight must be a positive integer')
        key = (ind, sig)
        if key in seen:
            raise ValueError('Duplicate signed atom is ambiguous')
        seen.add(key)
        q = F(w, denominator)
        for x, i in enumerate(ind):
            for y, j in enumerate(ind):
                r[i][j] -= q * sig[x] * sig[y]
    return r, metric


def verify_block(a, item, gamma,allow_amplitudes=False):
    r, metric = residual(a, item.get('metric_weights'), item.get('atoms'),
                         item.get('atom_scale'),allow_amplitudes)
    margins = [r[i][i] - sum(abs(r[i][j]) for j in range(len(r)) if j != i)
               for i in range(len(r))]
    bound = min(v / m for v, m in zip(margins, metric))
    target = F(gamma)
    if bound < target:
        raise ValueError('Signed atoms leave a residual below the proposed threshold')
    return {'dimension': len(r), 'atom_count': len(item['atoms']),
            'atom_maximum_support': max([0] + [len(x['indices']) for x in item['atoms']]),
            'packing_threshold_lower': rational_text(bound),
            'packing_threshold_lower_float': float(bound),
            'minimum_diagonal_margin': rational_text(min(m / q for m, q in zip(margins, metric)) - target)}


def replay(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    rational=certificate.get('kind')=='spin_rational_atom_complement_v1'
    if not rational and certificate.get('kind') != 'spin_signed_atom_complement_v1':
        raise ValueError('Unsupported signed atom certificate')
    oracle = SpinZeroOracle(certificate)
    blocks = certificate.get('blocks')
    matrices = checked_blocks(oracle, certificate.get('retained_states'), blocks,max_block_dimension=384)
    gamma = F(certificate['target_lower'])
    rows = [verify_block(a, b, gamma,rational) for a, b in zip(matrices, blocks)]
    receipt={'target_lower': rational_text(gamma), 'target_lower_float': float(gamma),
            'blocks': rows, 'unique_action_states': len(oracle.cache),
            'referenced_determinants': oracle.referenced_state_count(),
            'target_certified': True,
            'scope': ('Exact Q positivity from complete physical block coverage, '
                      'arbitrary signed equal-amplitude support-three/four rank-one atoms, '
                      'and an exactly diagonally-dominant residual in the declared metric. '
                      'This does not certify arbitrary amplitudes, dense PSD residuals, '
                      'other metrics, or general quantum marginal representability.')}
    if rational:
        receipt['scope']='Exact Q positivity from complete physical block coverage, explicit rational-amplitude rank-one atoms on at most four coordinates, and a complete exact diagonally dominant residual against gamma W². This certifies this finite decomposition, not sufficiency of all factor-width-four atoms or scalable general marginal representability.'
    return receipt


def price_atoms(b, known, seed=0, batch=256,vector_mode='signed',pair_pricing='sampled',diagnostics=None,batch_selection='lowest'):
    """Bounded deterministic sparse-vector pricing; no completeness claim."""
    import numpy as np
    if pair_pricing not in ('sampled','all_pairs','all_supports') or (pair_pricing!='sampled' and vector_mode!='rational'):
        raise ValueError('Complete pair pricing requires rational vectors')
    if batch_selection not in ('lowest','by_support') or type(batch) is not int or batch<1:
        raise ValueError('Positive batch size and known selection rule required')
    n=len(b);supports=set();rng=np.random.default_rng(seed);coverage=None
    sizes=(2,3,4) if vector_mode=='rational' else (3,4)
    if pair_pricing=='all_supports':
        from experiments.marginal_complete_pricing import scan_supports
        coverage=scan_supports(b,top_k=max(256,4*batch))
        supports.update(tuple(x['support']) for entries in coverage['negative_supports_by_order'].values() for x in entries)
    if pair_pricing=='all_pairs': supports.update(combinations(range(n),2))
    for i in range(n):
        order=sorted((j for j in range(n) if j!=i),key=lambda j:(-abs(b[i,j]),j))[:6]
        for k in sizes:
            for tail in combinations(order,k-1): supports.add(tuple(sorted((i,)+tail)))
    if n>=3:
        _,vectors=np.linalg.eigh(b)
        for vector in vectors[:,:min(3,n)].T:
            top=np.argsort(-np.abs(vector),kind='stable')[:min(8,n)]
            for k in sizes:
                for group in combinations(top,k): supports.add(tuple(sorted(int(i) for i in group)))
        for _ in range(256):
            for k in sizes:
                if n>=k: supports.add(tuple(sorted(int(i) for i in rng.choice(n,k,replace=False))))
    violations=[];minimum=min(0.,coverage['minimum_eigenvalue']) if coverage else 0.;by_support={}
    for k in sizes:
        selected=sorted(s for s in supports if len(s)==k)
        if not selected: continue
        indices=np.asarray(selected);blocks=b[indices[:,:,None],indices[:,None,:]]
        if vector_mode=='rational':
            eigenvalues,eigenvectors=np.linalg.eigh(blocks)
            minimum=min(minimum,float(eigenvalues[:,0].min()))
            by_support[str(k)]={'supports':len(selected),'minimum_eigenvalue':float(eigenvalues[:,0].min()),
                               'negative_eigenspaces':int(np.count_nonzero(eigenvalues[:,0] < -1e-8))}
            for q in np.flatnonzero(eigenvalues[:,0] < -1e-8):
                v=eigenvectors[q,:,0];integer=np.rint(v/np.max(np.abs(v))*10**6).astype(int)
                keep=np.flatnonzero(integer)
                if len(keep)<2: continue
                integer=integer[keep]
                if integer[0]<0: integer=-integer
                scale=int(np.max(np.abs(integer)));group=tuple(selected[int(q)][int(i)] for i in keep)
                vector=tuple(F(int(x),scale) for x in integer);key=(group,vector)
                vf=np.array(vector,dtype=float);pairing=float(vf@b[np.ix_(group,group)]@vf)
                if pairing < -1e-8 and key not in known: violations.append((pairing,key))
            continue
        for bits in range(1<<(k-1)):
            signs=(1,)+tuple(1 if bits&(1<<q) else -1 for q in range(k-1))
            v=np.asarray(signs);values=np.einsum('i,nij,j->n',v,blocks,v)
            minimum=min(minimum,float(values.min()))
            for q in np.flatnonzero(values < -1e-8):
                key=(selected[int(q)],signs)
                if key not in known: violations.append((float(values[q]),key))
    # Eigenvectors on different principal supports can round to the same
    # smaller support. Load each canonical direction only once.
    unique={}
    for value,key in violations: unique[key]=min(value,unique.get(key,value))
    violations=sorted(((value,key) for key,value in unique.items()),key=lambda x:(x[0],x[1]))
    chosen=violations[:batch]
    if batch_selection=='by_support':
        chosen=[];quota=batch//len(sizes)
        for k in sizes: chosen.extend([entry for entry in violations if len(entry[1][0])==k][:quota])
        selected_keys={key for _,key in chosen}
        chosen.extend([entry for entry in violations if entry[1] not in selected_keys][:batch-len(chosen)])
        chosen.sort(key=lambda x:(x[0],x[1]))
    if diagnostics is not None:
        diagnostics.update(pair_pricing=pair_pricing,batch_selection=batch_selection,by_support=by_support,
                           all_pair_supports_checked=pair_pricing in ('all_pairs','all_supports'),
                           selected_by_support={str(k):sum(len(key[0])==k for _,key in chosen) for k in sizes},
                           scope='Numerical local eigenvalue diagnostics. Complete pair support coverage does not certify exact positivity or complete support-three/four pricing.')
        if coverage:
            diagnostics['complete_support_scan']={k:v for k,v in coverage.items() if k not in ('negative_supports','negative_supports_by_order')}
            diagnostics['scope']='Every support of size two through four is scanned numerically, with near-boundary cases recorded; a bounded set of negative supports plus heuristic candidates proposes rational cuts. This is not exact dual positivity or global cone membership.'
    return [key for _,key in chosen],minimum,coverage['checked_total'] if coverage else len(supports)


def propose(a,weights,seconds=180,round_cap=12,seed=0,initial_item=None,pair_mode='all',vector_mode='signed',pair_pricing='sampled',initial_dictionary=None,batch_selection='lowest',initial_basis=None):
    """All-pair dual LP with adaptive signed-atom rows and exact primal export."""
    import math
    import highspy
    import numpy as np
    from scipy.sparse import coo_matrix
    from experiments.marginal_clique_gap import candidates,balanced_signs
    if not math.isfinite(seconds) or seconds<=0 or type(round_cap) is not int or not 1<=round_cap<=32:
        raise ValueError('Positive finite time and bounded round budget required')
    if pair_mode not in ('all','active'): raise ValueError('Pair mode must be all or active')
    if vector_mode not in ('signed','rational'): raise ValueError('Unknown local vector family')
    if pair_pricing not in ('sampled','all_pairs','all_supports') or (pair_pricing!='sampled' and vector_mode!='rational'):
        raise ValueError('Complete pair pricing requires rational vectors')
    if batch_selection not in ('lowest','by_support'): raise ValueError('Unknown batch selection rule')
    rational=vector_mode=='rational'
    started=time.monotonic();n=len(a);h,metric=metric_matrix(a,weights)
    hf=np.array(h,dtype=float);mf=np.array(metric,dtype=float)
    initial_bound=None
    if initial_item is None:
        atoms=[(tuple(g),tuple(balanced_signs(h,g))) for g in candidates(h,4,min_edge=F(0))]
    else:
        if initial_item.get('metric_weights')!=weights: raise ValueError('Initial atoms must use the same metric')
        initial_residual,initial_metric=residual(a,weights,initial_item.get('atoms'),initial_item.get('atom_scale'),rational)
        initial_bound=min((initial_residual[i][i]-sum(abs(initial_residual[i][j]) for j in range(n) if j!=i))/initial_metric[i] for i in range(n))
        atoms=[atom_vector(atom,n,rational) for atom in initial_item['atoms']]
    initial_proof=list(atoms)
    if initial_dictionary is not None:
        if type(initial_dictionary) is not list or len(initial_dictionary)>100000:
            raise ValueError('Bounded explicit direction dictionary required')
        dictionary=[]
        for entry in initial_dictionary:
            if type(entry) is not dict or 'weight' in entry:
                raise ValueError('Dictionary entries must be unweighted vectors')
            dictionary.append(atom_vector(entry,n,rational))
        if len(set(dictionary))!=len(dictionary): raise ValueError('Duplicate dictionary direction')
        if not set(initial_proof)<=set(dictionary): raise ValueError('Dictionary must retain every initial proof direction')
        atoms=dictionary
    initial_count=len(atoms)
    known=set(atoms)
    required={e for group,signs in atoms for e in combinations(group,2)}
    pairs=[(i,j) for i,j in combinations(range(n),2) if pair_mode=='all' or h[i][j] or (i,j) in required]
    initial_pair_count=len(pairs);index={e:n+k for k,e in enumerate(pairs)}
    orientation={e:1 if h[e[0]][e[1]]>=0 else -1 for e in pairs}
    def atom_coefficients(atom):
        group,signs=atom;node={i:v*v for i,v in zip(group,signs)};entry={}
        for q,i in enumerate(group):
            for r in range(q+1,len(group)):
                j=group[r];value=orientation[i,j]*signs[q]*signs[r]
                node[i]-=value;node[j]-=value;entry[index[i,j]]=float(value)
        entry.update({i:float(F(v)/metric[i]) for i,v in node.items() if v})
        return entry
    entries=[{i:1. for i in range(n)}]
    entries += [{i:float(F(2)/metric[i]),j:float(F(2)/metric[j]),index[i,j]:-1.} for i,j in pairs]
    entries += [atom_coefficients(atom) for atom in atoms]
    atom_rows=list(range(1+len(pairs),len(entries)))
    # Stable semantic names let a native basis survive insertion of cuts.
    def atom_key(atom):
        return 'atom:'+json.dumps(encode_vector(atom),sort_keys=True,separators=(',',':'))
    column_keys=['node:'+str(i) for i in range(n)]+['pair:%d,%d'%e for e in pairs]
    row_keys=['normalization']+['paircap:%d,%d'%e for e in pairs]+[atom_key(x) for x in atoms]
    rr=[];cc=[];vv=[]
    for row,entry in enumerate(entries):
        for col,value in entry.items(): rr.append(row);cc.append(col);vv.append(value)
    matrix=coo_matrix((vv,(rr,cc)),shape=(len(entries),n+len(pairs))).tocsc()
    baseline=[h[i][i]-sum(abs(h[i][j]) for j in range(n) if j!=i) for i in range(n)]
    lp=highspy.HighsLp();lp.num_col_=matrix.shape[1];lp.num_row_=matrix.shape[0]
    lp.col_cost_=np.array([v/m for v,m in zip(baseline,metric)]+[abs(h[i][j]) for i,j in pairs],dtype=float)
    lp.col_lower_=np.zeros(matrix.shape[1]);lp.col_upper_=np.full(matrix.shape[1],highspy.kHighsInf)
    lp.row_lower_=np.array([1.]+[0.]*(matrix.shape[0]-1));lp.row_upper_=np.array([1.]+[highspy.kHighsInf]*(matrix.shape[0]-1))
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    solver=highspy.Highs()
    for key,value in [('output_flag',False),('threads',1),('primal_feasibility_tolerance',1e-9),('dual_feasibility_tolerance',1e-9)]:
        if solver.setOptionValue(key,value)!=highspy.HighsStatus.kOk: raise ValueError('Native option rejected')
    if solver.passModel(lp)==highspy.HighsStatus.kError: raise ValueError('Native atom LP rejected')
    applied_basis=False
    remapped_basis=None
    if initial_basis is not None:
        if type(initial_basis) is not dict or initial_basis.get('version')!=1 or initial_basis.get('dimension')!=n or initial_basis.get('metric_weights')!=weights or initial_basis.get('vector_mode')!=vector_mode or initial_basis.get('pair_mode')!=pair_mode:
            raise ValueError('Native basis metadata does not match this model')
        oldc=initial_basis.get('column_keys'); oldr=initial_basis.get('row_keys')
        cs=initial_basis.get('column_status'); rs=initial_basis.get('row_status')
        if (type(oldc) is not list or type(oldr) is not list or type(cs) is not list or type(rs) is not list
            or any(type(k) is not str for k in oldc+oldr)
            or len(set(oldc))!=len(oldc) or len(set(oldr))!=len(oldr)
            or len(cs)!=len(oldc) or len(rs)!=len(oldr)
            or not set(oldc)<=set(column_keys) or not set(oldr)<=set(row_keys)
            or 'normalization' not in oldr or not {'node:'+str(i) for i in range(n)}<=set(oldc)):
            raise ValueError('Native basis keys do not match retained model')
        try:
            cmap=dict(zip(oldc,cs)); rmap=dict(zip(oldr,rs))
            def status(value):
                if type(value) is not int or value not in range(5): raise ValueError('Invalid native basis status')
                return highspy.HighsBasisStatus(value)
            cstatus=[status(cmap[k]) if k in cmap else highspy.HighsBasisStatus.kLower for k in column_keys]
            rstatus=[status(rmap[k]) if k in rmap else highspy.HighsBasisStatus.kBasic for k in row_keys]
            if sum(x==highspy.HighsBasisStatus.kBasic for x in cstatus+rstatus)!=len(row_keys):
                raise ValueError('Native basis must have one basic variable per row')
            remapped_basis=(cstatus,rstatus)
        except (TypeError,ValueError):
            raise ValueError('Malformed native basis statuses')
    initial=highspy.HighsSolution();initial.col_value=np.r_[mf/mf.sum(),np.full(len(pairs),2/mf.sum())];initial.value_valid=True
    solver.setSolution(initial)
    if remapped_basis is not None:
        basis=highspy.HighsBasis();basis.col_status,basis.row_status=remapped_basis
        if solver.setBasis(basis)==highspy.HighsStatus.kError: raise ValueError('Native basis rejected')
        applied_basis=True
    best=[];best_bound=float(min(v/m for v,m in zip(baseline,metric)));history=[];reason='round_budget';status=None
    if initial_item is not None:
        best_bound=float(initial_bound)
        best=[(atom,float(F(data['weight'],initial_item['atom_scale']))) for atom,data in zip(initial_proof,initial_item['atoms'])]
    last_priced=[]
    for iteration in range(round_cap):
        remaining=seconds-(time.monotonic()-started)
        if remaining<=0: reason='time_budget';break
        solver.setOptionValue('time_limit',solver.getRunTime()+remaining)
        prior_runtime=solver.getRunTime()
        solver.run();status=solver.getModelStatus();sol=solver.getSolution()
        numerical_bound=None
        if sol.dual_valid:
            coefficients=np.maximum(0.,np.asarray(sol.row_dual)[atom_rows])
            r=hf.copy();current=[]
            for atom,value in zip(atoms,coefficients):
                if value>0:
                    group,signs=atom;v=np.asarray(signs,dtype=float);r[np.ix_(group,group)]-=value*np.outer(v,v)
                    current.append((atom,float(value)))
            numerical_bound=float(np.min((np.diag(r)-(np.abs(r).sum(axis=1)-np.abs(np.diag(r))))/mf))
            if numerical_bound>best_bound: best_bound=numerical_bound;best=current
        row={'round':iteration,'selected_atoms':len(atoms),'active_atoms':len(best),
             'active_pairs':len(pairs),'lp_seconds':solver.getRunTime()-prior_runtime,
             'numerical_primal_threshold':numerical_bound,'best_numerical_threshold':best_bound,
             'solver_status':solver.modelStatusToString(status)}
        if status!=highspy.HighsModelStatus.kOptimal or not sol.value_valid:
            history.append(row);reason=solver.modelStatusToString(status);break
        x=np.asarray(sol.col_value);y=x[:n]/mf;b=np.diag(y)
        for e,col in index.items():
            i,j=e;b[i,j]=b[j,i]=orientation[e]*(x[col]-y[i]-y[j])/2
        pricing={}
        new,minimum,count=price_atoms(b,known,seed+iteration,vector_mode=vector_mode,pair_pricing=pair_pricing,diagnostics=pricing,batch_selection=batch_selection)
        last_priced=new
        row.update(minimum_heuristic_pairing=minimum,priced_supports=count,added_atoms=len(new),pricing=pricing)
        history.append(row);print(json.dumps(row),flush=True)
        if not new:
            reason='complete_scan_no_new_rounded_cut' if pair_pricing=='all_supports' else 'heuristic_pricing_exhausted'
            break
        if iteration+1==round_cap: break
        missing=sorted({e for group,signs in new for e in combinations(group,2)}-index.keys())
        if missing:
            if pair_mode!='active' or any(h[i][j] for i,j in missing):
                raise ValueError('Only unused zero-Hamiltonian pairs may be activated')
            loaded=solver.addCols(len(missing),np.zeros(len(missing)),np.zeros(len(missing)),
                np.full(len(missing),highspy.kHighsInf),0,np.zeros(len(missing)+1,dtype=np.int32),
                np.array([],dtype=np.int32),np.array([],dtype=float))
            if loaded==highspy.HighsStatus.kError: raise ValueError('Native new pair columns rejected')
            for e in missing:
                index[e]=n+len(pairs)
                pairs.append(e)
                orientation[e]=1
                column_keys.append('pair:%d,%d'%e)
                row_keys.append('paircap:%d,%d'%e)
        additions=[{i:float(F(2)/metric[i]),j:float(F(2)/metric[j]),index[i,j]:-1.} for i,j in missing]
        row_start=solver.getNumRow()
        additions += [atom_coefficients(atom) for atom in new];start=[0];cols=[];values=[]
        for entry in additions:
            for col,value in sorted(entry.items()): cols.append(col);values.append(value)
            start.append(len(cols))
        loaded=solver.addRows(len(additions),np.zeros(len(additions)),np.full(len(additions),highspy.kHighsInf),len(cols),np.asarray(start,dtype=np.int32),np.asarray(cols,dtype=np.int32),np.asarray(values))
        if loaded==highspy.HighsStatus.kError: raise ValueError('Native priced rows rejected')
        atom_rows.extend(range(row_start+len(missing),row_start+len(additions)))
        row_keys.extend(atom_key(x) for x in new)
        atoms.extend(new);known.update(new)
    scale=10**15
    exported=[]
    for (group,vector),value in best:
        weight=int(value*scale)
        if weight<=0: continue
        atom=dict(encode_vector((group,vector)),weight=weight)
        exported.append(atom)
    item={'metric_weights':weights,'atom_scale':scale,'atoms':exported}
    r,m=residual(a,weights,exported,scale,rational)
    bound=min((r[i][i]-sum(abs(r[i][j]) for j in range(n) if j!=i))/m[i] for i in range(n))
    if initial_bound is not None and bound<initial_bound:
        item={k:initial_item[k] for k in ('metric_weights','atom_scale','atoms')}
        bound=initial_bound
    verify_block(a,item,bound,rational)
    native_basis=None
    basis=solver.getBasis()
    if basis.valid:
        if len(basis.col_status)!=len(column_keys) or len(basis.row_status)!=len(row_keys):
            raise ValueError('Native basis dimensions do not match semantic keys')
        native_basis={'version':1,'dimension':n,'metric_weights':weights,'pair_mode':pair_mode,'vector_mode':vector_mode,
                      'column_keys':column_keys,'row_keys':row_keys,
                      'column_status':[int(x) for x in basis.col_status],
                      'row_status':[int(x) for x in basis.row_status]}
    return item,bound,{'dimension':n,'all_pairs':n*(n-1)//2,'active_pairs':len(pairs),
        'initial_pairs':initial_pair_count,'pair_mode':pair_mode,'vector_mode':vector_mode,'pair_pricing':pair_pricing,'batch_selection':batch_selection,'candidate_atoms':len(atoms),
        'initial_atoms':initial_count,'initial_proof_atoms':len(initial_proof) if initial_item else 0,
        'atom_dictionary':[encode_vector(key) for key in dict.fromkeys(atoms+last_priced)],
        'initial_exact_threshold':rational_text(initial_bound) if initial_bound is not None else None,
        'active_atoms':len(item['atoms']),'stopping_reason':reason,'history':history,
        'elapsed_seconds':time.monotonic()-started,'seconds_cap':seconds,'round_cap':round_cap,
        'exact_threshold':rational_text(bound),'exact_threshold_float':float(bound),'native_basis':native_basis,
        'initial_basis_applied':applied_basis,
        'scope':'Adaptive dual LP; every nonzero Hamiltonian pair and every pair used by a selected atom is modeled. Inactive pairs have identically zero primal residual; the numerical pricing dual is extended there by zero. Row multipliers propose positive atoms. Bounded support pricing is heuristic, not a global optimality test. The complete exact DD residual determines the exported threshold.'}


def construct(source,output,seconds=180,round_cap=12,atom_source=None,pair_mode='all',vector_mode='signed',pair_pricing='sampled',dictionary_source=None,batch_selection='lowest'):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    old=json.loads(Path(source).read_text());out=Path(output)
    if out.exists(): raise ValueError('Preserve previous signed atom export')
    if old.get('kind')!='spin_clique_metric_probe_v1': raise ValueError('Metric source required')
    base={k:old[k] for k in ('modes','particles','hamiltonian','retained_states')}
    initial=None
    if atom_source is not None:
        initial=json.loads(Path(atom_source).read_text())
        allowed=('spin_signed_atom_complement_v1','spin_rational_atom_complement_v1') if vector_mode=='rational' else ('spin_signed_atom_complement_v1',)
        if (initial.get('kind') not in allowed
                or any(initial.get(k)!=v for k,v in base.items())
                or [(b.get('states'),b.get('metric_weights')) for b in initial.get('blocks',[])]
                   !=[(b['states'],b['metric_weights']) for b in old['blocks']]):
            raise ValueError('Initial atoms must bind the same Hamiltonian, P, ordered blocks, and metrics')
    dictionary=None
    if dictionary_source is not None:
        dictionary=json.loads(Path(dictionary_source).read_text())
        if (type(dictionary) is not dict or type(dictionary.get('blocks')) is not list
                or any(type(b) is not dict or type(b.get('atom_dictionary')) is not list for b in dictionary['blocks'])):
            raise ValueError('Explicit block direction dictionaries required')
        if (dictionary.get('kind')!='spin_atom_dictionary_v1' or dictionary.get('vector_mode')!=vector_mode
                or any(dictionary.get(k)!=v for k,v in base.items())
                or [(b.get('states'),b.get('metric_weights')) for b in dictionary.get('blocks',[])]
                   !=[(b['states'],b['metric_weights']) for b in old['blocks']]):
            raise ValueError('Dictionary must bind the same vector family, Hamiltonian, P, ordered blocks, and metrics')
    oracle=SpinZeroOracle(base);matrices=checked_blocks(oracle,base['retained_states'],old['blocks'],max_block_dimension=384)
    blocks=[];bounds=[];history=[];dictionary_blocks=[]
    for i,(a,b) in enumerate(zip(matrices,old['blocks'])):
        item,bound,diagnostic=propose(a,b['metric_weights'],seconds,round_cap,initial_item=initial['blocks'][i] if initial else None,pair_mode=pair_mode,vector_mode=vector_mode,pair_pricing=pair_pricing,
            initial_dictionary=dictionary['blocks'][i]['atom_dictionary'] if dictionary else None,batch_selection=batch_selection,
            initial_basis=dictionary['blocks'][i].get('native_basis') if dictionary else None)
        dictionary_blocks.append({'states':b['states'],'metric_weights':b['metric_weights'],
                                  'atom_dictionary':diagnostic.pop('atom_dictionary'),
                                  'native_basis':diagnostic.pop('native_basis')})
        diagnostic['discovered_atoms']=len(dictionary_blocks[-1]['atom_dictionary'])
        blocks.append(dict(states=b['states'],**item));bounds.append(bound);history.append(diagnostic)
        print(json.dumps(diagnostic),flush=True)
    c=dict(base,kind='spin_rational_atom_complement_v1' if vector_mode=='rational' else 'spin_signed_atom_complement_v1',blocks=blocks,target_lower=rational_text(min(bounds)))
    receipt=replay(c)
    receipt.update(requested_target_lower=old['target_lower'],requested_target_certified=min(bounds)>=F(old['target_lower']),
                   initial_atom_source=str(atom_source) if atom_source is not None else None,
                   initial_dictionary_source=str(dictionary_source) if dictionary_source is not None else None)
    search_state=dict(base,kind='spin_atom_dictionary_v1',vector_mode=vector_mode,blocks=dictionary_blocks,
                      scope='Discovered direction dictionary plus an optional semantic native LP basis. Basis is numerical warm-start state only; exact replay remains the certificate.')
    out.mkdir(parents=True)
    for name,data in [('certificate.json',c),('receipt.json',receipt),('proposal.json',history),('search_state.json',search_state)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--verify');p.add_argument('--source');p.add_argument('--output')
    p.add_argument('--seconds',type=float,default=180);p.add_argument('--round-cap',type=int,default=12);p.add_argument('--atom-source');p.add_argument('--pair-mode',choices=['all','active'],default='all');p.add_argument('--vector-mode',choices=['signed','rational'],default='signed');p.add_argument('--pair-pricing','--support-pricing',choices=['sampled','all_pairs','all_supports'],default='sampled');p.add_argument('--dictionary-source');p.add_argument('--batch-selection',choices=['lowest','by_support'],default='lowest');args=p.parse_args()
    if args.verify: print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
    elif args.source and args.output: print(json.dumps(construct(args.source,args.output,args.seconds,args.round_cap,args.atom_source,args.pair_mode,args.vector_mode,args.pair_pricing,args.dictionary_source,args.batch_selection),indent=2))
    else: p.error('Specify --verify or --source and --output')
