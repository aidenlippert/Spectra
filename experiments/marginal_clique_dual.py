"""Exact upper bounds on the fixed-metric balanced-clique packing family.

This separates a restricted square family, not the physical PSD cone or FW4.
Numerical proposals are rounded and repaired; replay uses only integer/rational
inequalities and enumerates every balanced clique of support three or four.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time

from experiments.marginal_clique_gap import metric_matrix,candidates
from experiments.marginal_implicit_certificate import rational_text


def geometry(a,weights):
    a,m=metric_matrix(a,weights);n=len(a)
    edges={(i,j):abs(a[i][j]) for i in range(n) for j in range(i+1,n) if a[i][j]}
    b=[a[i][i]-sum(abs(a[i][j]) for j in range(n) if j!=i) for i in range(n)]
    return a,m,edges,b,candidates(a,4,min_edge=F(0))


def verify_block(a,item):
    a,m,edges,b,groups=geometry(a,item.get('metric_weights'))
    y,z=item.get('node_weights'),item.get('edge_weights')
    for values,length in ((y,len(a)),(z,len(edges))):
        if (type(values) is not list or len(values)!=length
                or any(type(v) is not int or not 0<=v<=10**24 for v in values)):
            raise ValueError('Bounded nonnegative integer dual weights required')
    normalization=sum(mi*yi for mi,yi in zip(m,y))
    if normalization<=0: raise ValueError('Nonzero node normalization required')
    lookup=dict(zip(edges,z))
    for group in groups:
        if sum(lookup[e] for e in combinations(group,2)) < (len(group)-2)*sum(y[i] for i in group):
            raise ValueError('A balanced clique violates the exact dual inequality')
    upper=(sum(x*yi for x,yi in zip(b,y))+sum(c*zi for c,zi in zip(edges.values(),z)))/normalization
    return {'dimension':len(a),'cliques_checked':len(groups),'edge_weights':len(z),
            'packing_threshold_upper':rational_text(upper),'packing_threshold_upper_float':float(upper)}


def replay(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    if certificate.get('kind')!='spin_clique_packing_dual_v1':
        raise ValueError('Unsupported clique packing dual')
    oracle=SpinZeroOracle(certificate);blocks=certificate.get('blocks')
    matrices=checked_blocks(oracle,certificate.get('retained_states'),blocks)
    rows=[verify_block(a,b) for a,b in zip(matrices,blocks)]
    upper=min(F(r['packing_threshold_upper']) for r in rows)
    target=F(certificate['target_lower'])
    return {'packing_threshold_upper':rational_text(upper),'packing_threshold_upper_float':float(upper),
        'target_lower':rational_text(target),'target_excluded':upper<target,'blocks':rows,
        'unique_action_states':len(oracle.cache),
        'scope':'Upper bound on the best complement threshold obtainable by balanced support-three/four clique squares plus sign-matched leftover edge squares in the declared fixed diagonal metrics. Every nonzero-edge balanced clique is checked, including cliques omitted from numerical pricing. Does not bound the physical spectrum from above or exclude general factor-width-four, cancellations, or alternative metrics.'}


def replay_star(certificate):
    """One-row dual: each incident triangle needs one cover edge, each quad two."""
    from experiments.marginal_spin_reduction import SpinZeroOracle
    if certificate.get('kind')!='spin_clique_star_obstruction_v1':
        raise ValueError('Unsupported local clique obstruction')
    oracle=SpinZeroOracle(certificate);p=oracle.retained(certificate.get('retained_states'))
    row=certificate.get('row_state')
    if not oracle.valid_state(row) or row in p: raise ValueError('Physical Q row required')
    action=oracle.action(row)
    states=certificate.get('metric_states')
    expected={row}|{s for s,v in action.items() if v and s not in p}
    if (type(states) is not list or not 1<=len(states)<=256 or states!=sorted(expected)):
        raise ValueError('Metric must cover exactly the complete Q row neighborhood')
    a=[[oracle.action(s).get(t,F(0)) for t in states] for s in states]
    if certificate.get('metric_weights') is None: raise ValueError('Explicit local metric required')
    a,m=metric_matrix(a,certificate['metric_weights'])
    index={s:i for i,s in enumerate(states)};i=index[row]
    cover=certificate.get('covered_edges')
    if type(cover) is not list or len(cover)>32640: raise ValueError('Bounded edge cover required')
    edges=set()
    for edge in cover:
        if (type(edge) is not list or len(edge)!=2
                or any(type(s) is not int or s not in index for s in edge)
                or edge!=sorted(set(edge)) or tuple(edge) in edges
                or not a[index[edge[0]]][index[edge[1]]]):
            raise ValueError('Distinct nonzero ordered neighborhood edges required')
        edges.add(tuple(edge))
    neighbors=[s for s in states if s!=row]
    signs={s:1 if a[i][index[s]]>0 else -1 for s in neighbors}
    def balanced(s,t):
        v=a[index[s]][index[t]]
        return bool(v) and (1 if v>0 else -1)==signs[s]*signs[t]
    def covered(group):
        return sum(tuple(sorted(e)) in edges for e in combinations(group,2))
    triangles=quads=0
    for s,t in combinations(neighbors,2):
        if balanced(s,t):
            triangles+=1
            if covered([row,s,t])<1: raise ValueError('Uncovered balanced triangle')
    for s,t,u in combinations(neighbors,3):
        if balanced(s,t) and balanced(s,u) and balanced(t,u):
            quads+=1
            if covered([row,s,t,u])<2: raise ValueError('Balanced quad needs two cover edges')
    upper=(a[i][i]-sum(abs(a[i][j]) for j in range(len(a)) if j!=i)
           +sum(abs(a[index[s]][index[t]]) for s,t in edges))/m[i]
    target=F(certificate['target_lower'])
    if upper>=target: raise ValueError('Local cover does not exclude the target')
    return {'packing_threshold_upper':rational_text(upper),'packing_threshold_upper_float':float(upper),
        'target_lower':rational_text(target),'target_excluded':True,'row_state':row,
        'metric_states':len(states),'covered_edges':len(edges),
        'balanced_triangles_checked':triangles,'balanced_quads_checked':quads,
        'unique_action_states':len(oracle.cache),'referenced_determinants':oracle.referenced_state_count(),
        'scope':'Exact one-row dual obstruction to balanced support-three/four clique packing plus sign-matched leftover edge squares. Every clique containing the active row is checked; all other clique inequalities have zero left side. Only the complete Q row neighborhood is enumerated. Applies to the declared metric weights and any positive extension outside that neighborhood. Does not exclude arbitrary local amplitudes, cancellation between differently signed squares, or other metrics.'}


def compress_star(source,output):
    old=json.loads(Path(source).read_text());out=Path(output)
    if out.exists(): raise ValueError('Preserve previous star export')
    # Compression is merely a proposal; standalone replay checks the original Hamiltonian.
    if old.get('kind')!='spin_clique_packing_dual_v1': raise ValueError('Packing dual required')
    from experiments.marginal_spin_reduction import SpinZeroOracle
    oracle=SpinZeroOracle(old);p=oracle.retained(old['retained_states'])
    for block in old['blocks']:
        y=block['node_weights'];nonzero=[i for i,v in enumerate(y) if v]
        if len(nonzero)!=1: continue
        active=nonzero[0];row=block['states'][active]
        states=block['states'];a=[[oracle.action(s).get(t,F(0)) for t in states] for s in states]
        pairs=[(s,t) for i,s in enumerate(states) for j,t in enumerate(states) if i<j and a[i][j]]
        if any(z not in (0,y[active]) for z in block['edge_weights']): continue
        cover=[sorted(e) for e,z in zip(pairs,block['edge_weights']) if z]
        local=sorted({row}|{s for s,v in oracle.action(row).items() if v and s not in p})
        lookup=dict(zip(states,block['metric_weights']))
        c={k:old[k] for k in ('modes','particles','hamiltonian','retained_states','target_lower')}
        c.update(kind='spin_clique_star_obstruction_v1',row_state=row,metric_states=local,
                 metric_weights=[lookup[s] for s in local],covered_edges=cover)
        receipt=replay_star(c);out.mkdir(parents=True)
        for name,data in [('certificate.json',c),('receipt.json',receipt)]:
            (out/name).write_text(json.dumps(data,indent=2)+'\n')
        return receipt
    raise ValueError('No one-row unit edge-cover dual available')


def propose(a,weights,seconds=60):
    import highspy
    import numpy as np
    from scipy.sparse import coo_matrix
    original=a
    a,m,edges,b,groups=geometry(a,weights);n=len(a)
    index={e:n+i for i,e in enumerate(edges)}
    rows=[0]*n;cols=list(range(n));values=[1.]*n
    for row,group in enumerate(groups,1):
        for i in group:
            rows.append(row);cols.append(i);values.append(float(F(2-len(group))/m[i]))
        for e in combinations(group,2):
            rows.append(row);cols.append(index[e]);values.append(1.)
    matrix=coo_matrix((values,(rows,cols)),shape=(len(groups)+1,n+len(edges))).tocsc()
    h=highspy.Highs()
    for k,v in [('output_flag',False),('threads',1),('time_limit',seconds),
                ('primal_feasibility_tolerance',1e-9),('dual_feasibility_tolerance',1e-9)]:
        if h.setOptionValue(k,v)!=highspy.HighsStatus.kOk: raise ValueError('LP option rejected')
    lp=highspy.HighsLp();lp.num_col_=matrix.shape[1];lp.num_row_=matrix.shape[0]
    lp.col_cost_=np.array([bi/mi for bi,mi in zip(b,m)]+list(edges.values()),dtype=float)
    lp.col_lower_=np.zeros(matrix.shape[1]);lp.col_upper_=np.full(matrix.shape[1],highspy.kHighsInf)
    lp.row_lower_=np.array([1.]+[0.]*len(groups));lp.row_upper_=np.array([1.]+[highspy.kHighsInf]*len(groups))
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    if h.passModel(lp)==highspy.HighsStatus.kError: raise ValueError('Dual LP rejected')
    # Known feasible star cover, also available if the numerical solver returns no proposal.
    best=min(range(n),key=lambda i:(b[i]+F(2,3)*sum(c for e,c in edges.items() if i in e))/m[i])
    fallback=np.zeros(matrix.shape[1]);fallback[best]=1.
    for e,col in index.items():
        if best in e: fallback[col]=float(F(2,3)/m[best])
    initial=highspy.HighsSolution();initial.col_value=fallback;initial.value_valid=True
    h.setSolution(initial)
    started=time.monotonic();h.run();solution=h.getSolution()
    x=np.array(solution.col_value) if solution.value_valid else fallback
    scale=10**9
    y=[max(0,round(float(x[i])/float(m[i])*scale)) for i in range(n)]
    if not any(y): y[best]=scale
    z={e:max(0,round(float(x[col])*scale)) for e,col in index.items()}
    repaired=0
    for group in groups:
        ce=list(combinations(group,2))
        deficit=(len(group)-2)*sum(y[i] for i in group)-sum(z[e] for e in ce)
        if deficit>0:
            edge=min(ce,key=lambda e:(edges[e],e))
            z[edge]+=deficit;repaired+=deficit
    recipe={'metric_weights':weights,'node_weights':y,'edge_weights':list(z.values())}
    receipt=verify_block(original,recipe)
    return recipe,{'solver_status':h.modelStatusToString(h.getModelStatus()),
        'solver_seconds':time.monotonic()-started,'numerical_objective':float(lp.col_cost_@x),
        'integer_repair_mass':repaired,**receipt}


def construct(source,output,seconds=60):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    old=json.loads(Path(source).read_text());out=Path(output)
    if out.exists(): raise ValueError('Preserve previous dual export')
    if old.get('kind')!='spin_clique_metric_probe_v1': raise ValueError('Metric probe required')
    base={k:old[k] for k in ('modes','particles','hamiltonian','retained_states','target_lower')}
    oracle=SpinZeroOracle(base);matrices=checked_blocks(oracle,base['retained_states'],old['blocks'])
    blocks=[];history=[]
    for a,b in zip(matrices,old['blocks']):
        recipe,diagnostic=propose(a,b['metric_weights'],seconds)
        blocks.append(dict(states=b['states'],**recipe));history.append(diagnostic)
        print(json.dumps(diagnostic),flush=True)
    certificate=dict(base,kind='spin_clique_packing_dual_v1',blocks=blocks)
    receipt=replay(certificate);out.mkdir(parents=True)
    for name,data in [('certificate.json',certificate),('receipt.json',receipt),('proposal.json',history)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify');p.add_argument('--source')
    p.add_argument('--verify-star');p.add_argument('--compress-star')
    p.add_argument('--output');p.add_argument('--seconds',type=float,default=60);args=p.parse_args()
    if args.verify: print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
    elif args.verify_star: print(json.dumps(replay_star(json.loads(Path(args.verify_star).read_text())),indent=2))
    elif args.compress_star and args.output: print(json.dumps(compress_star(args.compress_star,args.output),indent=2))
    elif args.source and args.output: construct(args.source,args.output,args.seconds)
    else: p.error('Specify --verify or --source and --output')
