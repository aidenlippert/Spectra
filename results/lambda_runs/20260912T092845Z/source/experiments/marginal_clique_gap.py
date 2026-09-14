"""Reference gaps from sign-consistent clique squares and leftover edge squares.

This is a restricted factor-width cone, not every width-four decomposition.
LP values only propose weights; exact edge budgets and diagonal margins decide.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time

from experiments.marginal_implicit_certificate import rational_text


def balanced_signs(a, indices):
    signs = [1]+[(a[indices[0]][j]>0)-(a[indices[0]][j]<0) for j in indices[1:]]
    if any(not x for x in signs):
        return None
    for i,j in combinations(range(len(indices)),2):
        value=a[indices[i]][indices[j]]
        if not value or ((value>0)-(value<0)) != signs[i]*signs[j]:
            return None
    return signs


def margins(a, atoms, denominator):
    n=len(a)
    if type(denominator) is not int or not 1<=denominator<=10**18:
        raise ValueError('Bounded positive integer clique scale required')
    if type(atoms) is not list or len(atoms)>100000:
        raise ValueError('Bounded explicit clique atoms required')
    remaining={(i,j):abs(a[i][j]) for i in range(n) for j in range(i+1,n) if a[i][j]}
    diagonal=[a[i][i]-sum(abs(a[i][j]) for j in range(n) if j!=i) for i in range(n)]
    supports=set()
    for atom in atoms:
        if type(atom) is not dict:
            raise ValueError('Explicit clique atom required')
        indices,weight=atom.get('indices'),atom.get('weight')
        if (type(indices) is not list or not 3<=len(indices)<=4
                or any(type(i) is not int or not 0<=i<n for i in indices)
                or indices!=sorted(set(indices)) or tuple(indices) in supports
                or type(weight) is not int or weight<=0):
            raise ValueError('Distinct sorted three/four-coordinate atoms and positive integer weights required')
        if balanced_signs(a,indices) is None:
            raise ValueError('Clique square signs do not match all Hamiltonian edges')
        supports.add(tuple(indices))
        value=F(weight,denominator)
        for i,j in combinations(indices,2):
            remaining[i,j]-=value
            if remaining[i,j]<0:
                raise ValueError('Clique allocation overspends an exact edge budget')
        for i in indices:
            diagonal[i]+=(len(indices)-2)*value
    return diagonal,remaining


def metric_matrix(a,weights=None):
    """Congruence WHW and its identity metric W², using exact rational W."""
    n=len(a)
    if weights is None:
        return a,[F(1)]*n
    if (type(weights) is not list or len(weights)!=n
            or any(type(x) is not int or not 0<x<=10**12 for x in weights)):
        raise ValueError('Bounded positive integer diagonal metric required')
    w=[F(x,max(weights)) for x in weights]
    return [[a[i][j]*w[i]*w[j] for j in range(n)] for i in range(n)],[x*x for x in w]


def verify_cliques(a,gamma,item):
    transformed,metric=metric_matrix(a,item.get('metric_weights'))
    diagonal,remaining=margins(transformed,item.get('clique_atoms'),item.get('clique_scale'))
    margin=min(x/m for x,m in zip(diagonal,metric))-F(gamma)
    if margin<0:
        raise ValueError('Clique/edge squares do not certify the proposed threshold')
    atoms=item['clique_atoms']
    return {'dimension':len(a),'margin':margin,'clique_count':len(atoms),
            'maximum_atom_support':max([2]+[len(x['indices']) for x in atoms]),
            'remaining_edge_squares':sum(x>0 for x in remaining.values())}


def replay(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    if certificate.get('kind')!='spin_clique_complement_v1':
        raise ValueError('Unsupported clique complement certificate')
    oracle=SpinZeroOracle(certificate)
    blocks=certificate.get('blocks')
    matrices=checked_blocks(oracle,certificate.get('retained_states'),blocks)
    gamma=F(certificate['complement_lower'])
    receipts=[verify_cliques(a,gamma,b) for a,b in zip(matrices,blocks)]
    return {'complement_lower':rational_text(gamma),'complement_lower_float':float(gamma),
        'blocks':[{k:rational_text(v) if isinstance(v,F) else v for k,v in r.items()} for r in receipts],
        'unique_action_states':len(oracle.cache),'referenced_determinants':oracle.referenced_state_count(),
        'scope':'Exact Q positivity from complete block coverage, balanced clique rank-one squares, nonnegative leftover edge squares, and nonnegative diagonal residuals. All Q configurations remain explicit. Atoms have support at most four in determinant coordinates; no dense positivity factor is needed in replay. Restricted packing cone, not full factor-width-four membership.'}


def candidates(a,width,min_edge=F(1,10**8)):
    if width not in (3,4):
        raise ValueError('Clique width must be three or four')
    n=len(a);adj=[{j for j in range(i+1,n) if a[i][j] and abs(a[i][j])>=min_edge} for i in range(n)]
    result=[]
    for i in range(n):
        for j in sorted(adj[i]):
            common=adj[i]&adj[j]
            for k in sorted(common):
                tri=[i,j,k]
                if balanced_signs(a,tri) is None:
                    continue
                result.append(tri)
                if width==4:
                    for ell in sorted(common&adj[k]):
                        group=tri+[ell]
                        if balanced_signs(a,group) is not None:
                            result.append(group)
                if len(result)>100000:
                    raise ValueError('Clique candidate budget exceeded')
    return result


def incremental_packing(matrix,rhs,seconds):
    """Retain a native LP basis while adding dual-priced candidate columns."""
    import highspy
    import numpy as np
    h=highspy.Highs()
    for key,value in [('output_flag',False),('threads',1),('primal_feasibility_tolerance',1e-10),
                      ('dual_feasibility_tolerance',1e-9),('small_matrix_value',1e-12)]:
        if h.setOptionValue(key,value)!=highspy.HighsStatus.kOk:
            raise ValueError('Rejected native LP option '+key)
    first=matrix[:,:1].tocsc();count=matrix.shape[1]
    lp=highspy.HighsLp();lp.num_col_=1;lp.num_row_=matrix.shape[0]
    lp.col_cost_=[-1.];lp.col_lower_=[0.];lp.col_upper_=[highspy.kHighsInf]
    lp.row_lower_=[-highspy.kHighsInf]*matrix.shape[0];lp.row_upper_=rhs
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_=first.indptr;lp.a_matrix_.index_=first.indices;lp.a_matrix_.value_=first.data
    if h.passModel(lp)==highspy.HighsStatus.kError: raise ValueError('Native LP model rejected')
    active=[0];selected=np.zeros(count,dtype=bool);selected[0]=True
    solution=np.zeros(count);dual=None;history=[];started=time.monotonic();reason='round_budget'
    status=None
    for iteration in range(33):
        h.setOptionValue('time_limit',h.getRunTime()+max(.01,seconds-(time.monotonic()-started)))
        h.run();status=h.getModelStatus();current=h.getSolution()
        if current.value_valid:
            solution[active]=np.asarray(current.col_value)
        if status!=highspy.HighsModelStatus.kOptimal:
            reason=h.modelStatusToString(status);break
        dual=np.asarray(current.row_dual)
        reduced=np.asarray(-matrix.T@dual).ravel();reduced[selected]=np.inf
        best=float(np.min(reduced)) if np.any(~selected) else 0.
        history.append({'round':iteration,'selected_cliques':len(active)-1,
                        'numerical_gain':float(solution[0]),'minimum_unselected_reduced_cost':best})
        if best>=-1e-8:
            reason='priced_dictionary_optimal';break
        if iteration==32 or time.monotonic()-started>=seconds:
            break
        eligible=np.flatnonzero(reduced < -1e-8)
        order=eligible[np.lexsort((eligible,reduced[eligible]))][:512]
        columns=matrix[:,order].tocsc()
        loaded=h.addCols(len(order),np.zeros(len(order)),np.zeros(len(order)),
            np.full(len(order),highspy.kHighsInf),len(columns.data),columns.indptr,columns.indices,columns.data)
        if loaded==highspy.HighsStatus.kError: raise ValueError('Native LP rejected priced columns')
        active.extend(int(i) for i in order);selected[order]=True
    return solution,dual,{'solver_status':str(status),'stopping_reason':reason,
                          'selected_cliques':len(active)-1,'lp_history':history,'lp_seconds':time.monotonic()-started}


def propose_block(a,width,seconds=180,metric_weights=None):
    import numpy as np
    from scipy.sparse import coo_matrix
    a,metric=metric_matrix(a,metric_weights)
    n=len(a);groups=candidates(a,width)
    edges={(i,j):abs(a[i][j]) for i in range(n) for j in range(i+1,n) if a[i][j]}
    edge_index={edge:n+k for k,edge in enumerate(edges)}
    baseline=[a[i][i]-sum(abs(a[i][j]) for j in range(n) if j!=i) for i in range(n)]
    origin=min(x/m for x,m in zip(baseline,metric))
    rows=list(range(n));cols=[0]*n;values=[1.]*n
    for column,group in enumerate(groups,1):
        for i in group: rows.append(i);cols.append(column);values.append(float(F(2-len(group))/metric[i]))
        for edge in combinations(group,2): rows.append(edge_index[edge]);cols.append(column);values.append(1.)
    constraints=coo_matrix((values,(rows,cols)),shape=(n+len(edges),len(groups)+1)).tocsc()
    print(json.dumps({'phase':'packing_lp','dimension':n,'candidates':len(groups),'edges':len(edges),'initial_threshold':float(origin)}),flush=True)
    solution,dual,solver=incremental_packing(constraints,np.array([x/m-origin for x,m in zip(baseline,metric)]+list(edges.values()),dtype=float),seconds)
    denominator=10**18 if metric_weights is not None else 10**12
    weights=[max(0,int(np.floor(x*denominator))) for x in solution[1:]]
    loads={edge:0 for edge in edges}
    incidence={edge:[] for edge in edges}
    for index,(group,w) in enumerate(zip(groups,weights)):
        if w:
            for edge in combinations(group,2): loads[edge]+=w;incidence[edge].append(index)
    removed=0
    for edge,capacity in edges.items():
        debt=loads[edge]-int(capacity*denominator)
        if debt<=0: continue
        for index in sorted(incidence[edge],key=lambda i:(-weights[i],i)):
            amount=min(weights[index],debt)
            weights[index]-=amount;debt-=amount;removed+=amount
            for other in combinations(groups[index],2): loads[other]-=amount
            if not debt: break
        if debt: raise ValueError('Exact edge allocation repair failed')
    atoms=[{'indices':group,'weight':w} for group,w in zip(groups,weights) if w]
    diagonal,_=margins(a,atoms,denominator)
    bound=min(x/m for x,m in zip(diagonal,metric))
    recipe={'clique_scale':denominator,'clique_atoms':atoms}
    if metric_weights is not None: recipe['metric_weights']=metric_weights
    return recipe,bound,{
        'dimension':n,'candidate_cliques':len(groups),'active_cliques':len(atoms),
        **solver,'numerical_threshold':float(origin)+float(solution[0]),'exact_export_threshold':rational_text(bound),
        'exact_export_threshold_float':float(bound),'removed_coefficient_mass':rational_text(F(removed,denominator))}


def construct(source,output,width=4,metric_source=None,seconds=180):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    source,out=Path(source),Path(output)
    if out.exists(): raise ValueError('Preserve previous clique export')
    old=json.loads(source.read_text())['spin_symmetric_certificate']
    base={key:old[key] for key in ('modes','particles','hamiltonian','retained_states')}
    oracle=SpinZeroOracle(base)
    # Import only the declared component states, never the old LDL entries.
    blocks=[{'states':b['states']} for b in old['blocks']]
    if metric_source is not None:
        metric=json.loads(Path(metric_source).read_text())
        if (metric.get('kind')!='spin_clique_metric_probe_v1'
                or any(metric.get(k)!=v for k,v in base.items())
                or [b.get('states') for b in metric.get('blocks',[])]!=[b['states'] for b in blocks]):
            raise ValueError('Metric must bind the same Hamiltonian, particle sector, P, and ordered blocks')
        for b,m in zip(blocks,metric['blocks']):
            if m.get('metric_weights') is None: raise ValueError('Explicit metric weights required')
            b['metric_weights']=m['metric_weights']
    matrices=checked_blocks(oracle,base['retained_states'],blocks)
    started=time.monotonic();history=[];bounds=[]
    for block,a in zip(blocks,matrices):
        recipe,bound,diagnostic=propose_block(a,width,seconds,block.get('metric_weights'))
        block.update(recipe);bounds.append(bound);history.append(diagnostic)
        print(json.dumps(diagnostic),flush=True)
    certificate=dict(base,kind='spin_clique_complement_v1',blocks=blocks,complement_lower=rational_text(min(bounds)))
    receipt=replay(certificate)
    receipt.update(source=str(source),elapsed_seconds=time.monotonic()-started,
        discovery_scope='Imports reference Hamiltonian, P, and component state lists only. Clique allocation uses a numerical LP; exact rational edge budgets and margins accept the result. Full configuration/block construction remains explicit.')
    out.mkdir(parents=True)
    for name,data in [('certificate.json',certificate),('receipt.json',receipt),('proposal.json',history)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'threshold':receipt['complement_lower_float'],'seconds':receipt['elapsed_seconds']}),flush=True)


def replay_flat(certificate):
    """A single row separates the cone of equal-magnitude sparse squares."""
    from experiments.marginal_spin_reduction import SpinZeroOracle
    if certificate.get('kind')!='spin_flat_factor_obstruction_v1':
        raise ValueError('Unsupported equal-magnitude factor obstruction')
    width=certificate.get('maximum_support')
    if type(width) is not int or not 2<=width<=64:
        raise ValueError('Bounded integer support limit required')
    oracle=SpinZeroOracle(certificate)
    retained=oracle.retained(certificate.get('retained_states'))
    state=certificate.get('row_state')
    if not oracle.valid_state(state) or state in retained:
        raise ValueError('A physical complementary row is required')
    row=oracle.action(state)
    diagonal=row.get(state,F(0))
    off=sum(abs(x) for s,x in row.items() if s!=state and s not in retained)
    ceiling=diagonal-off/(width-1)
    target=F(certificate['target_lower'])
    if target<=ceiling:
        raise ValueError('This row does not obstruct the equal-magnitude cone at the target')
    return {'maximum_support_excluded':width,'target_lower':rational_text(target),
        'threshold_ceiling':rational_text(ceiling),'threshold_ceiling_float':float(ceiling),
        'normalized_dual_pairing':rational_text(ceiling-target),
        'row_diagonal':rational_text(diagonal),'q_offdiagonal_absolute_sum':rational_text(off),
        'unique_action_states':len(oracle.cache),
        'witness_support':1+sum(s!=state and s not in retained for s in row),
        'scope':'Any sum of PSD rank-one atoms whose nonzero amplitudes have equal magnitude and support at most k satisfies A_ii >= sum_j!=i |A_ij|/(k-1). The exact row violates this necessary condition for QHQ-target. Includes arbitrary signs and coefficient cancellations. Does not exclude arbitrary-amplitude factor-width-k decompositions.'}


def row_ceiling(a,weights,width=4):
    n=len(a)
    if (width not in (3,4) or type(weights) is not list or len(weights)!=n
            or any(type(x) is not int or not 0<x<=10**12 for x in weights)):
        raise ValueError('Bounded positive integer diagonal metric required')
    return min(a[i][i]-sum(abs(a[i][j])*weights[j] for j in range(n) if j!=i)/((width-1)*weights[i]) for i in range(n))


def replay_metric(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    if certificate.get('kind')!='spin_clique_metric_probe_v1':
        raise ValueError('Unsupported clique metric probe')
    oracle=SpinZeroOracle(certificate);blocks=certificate.get('blocks')
    matrices=checked_blocks(oracle,certificate.get('retained_states'),blocks)
    target=F(certificate['target_lower']);rows=[]
    for a,b in zip(matrices,blocks):
        unweighted=row_ceiling(a,[1]*len(a))
        weighted=row_ceiling(a,b.get('metric_weights'))
        rows.append({'dimension':len(a),'uniform_threshold_ceiling':rational_text(unweighted),
            'rescaled_threshold_ceiling':rational_text(weighted),'rescaled_ceiling_float':float(weighted),
            'necessary_row_gate_passes':weighted>=target})
    return {'target_lower':rational_text(target),'blocks':rows,'unique_action_states':len(oracle.cache),
        'scope':'Exact evaluation of a necessary row condition after positive diagonal congruence W. Flat atoms in WHW correspond to reciprocal-weight unequal-amplitude atoms in H. Passing this condition proves neither cone membership nor a positive complement gap. All reference configurations are still represented.'}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify');p.add_argument('--verify-flat');p.add_argument('--verify-metric');p.add_argument('--source');p.add_argument('--output');p.add_argument('--width',type=int,default=4)
    p.add_argument('--metric-source');p.add_argument('--seconds',type=float,default=180)
    args=p.parse_args()
    if args.verify: print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
    elif args.verify_flat: print(json.dumps(replay_flat(json.loads(Path(args.verify_flat).read_text())),indent=2))
    elif args.verify_metric: print(json.dumps(replay_metric(json.loads(Path(args.verify_metric).read_text())),indent=2))
    elif args.source and args.output: construct(args.source,args.output,args.width,args.metric_source,args.seconds)
    else:p.error('Specify --verify or --source and --output')
