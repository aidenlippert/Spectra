"""Exact fixed-metric clique duals with unrestricted signs of residual edges.

Balanced three/four-coordinate squares may overspend Hamiltonian edges, with
opposite-sign edge squares compensating. This is still a restricted atom family.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from itertools import combinations

from experiments.marginal_clique_dual import geometry,verify_block
from experiments.marginal_implicit_certificate import rational_text


def verify_signed(a,item):
    receipt=verify_block(a,item)
    transformed,m,edges,b,groups=geometry(a,item.get('metric_weights'))
    y,z=item['node_weights'],item['edge_weights']
    for (i,j),value in zip(edges,z):
        if value>2*(y[i]+y[j]):
            raise ValueError('Opposite-sign edge square violates the dual')
    receipt['both_sign_edges_checked']=len(edges)
    return receipt


def dual_matrix(a,item):
    """Unnormalized dual in the WHW chart; zero on absent Hamiltonian edges."""
    transformed,m,edges,b,groups=geometry(a,item.get('metric_weights'))
    n=len(a);y,z=item['node_weights'],item['edge_weights']
    matrix=[[F(0) for j in range(n)] for i in range(n)]
    for i in range(n): matrix[i][i]=F(y[i])
    for (i,j),value in zip(edges,z):
        sign=1 if transformed[i][j]>0 else -1
        matrix[i][j]=matrix[j][i]=F(sign*(value-y[i]-y[j]),2)
    return matrix,m


def separator_pairing(a,item,indices,signs):
    if (type(indices) is not list or not 3<=len(indices)<=4
            or any(type(i) is not int or not 0<=i<len(a) for i in indices)
            or indices!=sorted(set(indices)) or type(signs) is not list
            or len(signs)!=len(indices) or any(type(s) is not int or s not in (-1,1) for s in signs)
            or signs[0]!=1):
        raise ValueError('Distinct sorted signed three/four-coordinate separator required')
    matrix,m=dual_matrix(a,item)
    value=sum(signs[i]*signs[j]*matrix[s][t]
              for i,s in enumerate(indices) for j,t in enumerate(indices))
    if value>=0: raise ValueError('Signed square does not violate the parent dual')
    norm=sum(F(1)/m[s] for s in indices)
    trace=sum(m[i]*matrix[i][i] for i in range(len(a)))
    return {'atom_support':len(indices),'unnormalized_dual_pairing':rational_text(value),
            'normalized_dual_atom_pairing':rational_text(value/(norm*trace)),
            'normalized_dual_atom_pairing_float':float(value/(norm*trace))}


def replay_separator(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    if certificate.get('kind')!='spin_signed_clique_separator_v1':
        raise ValueError('Unsupported signed-clique separator')
    parent=certificate.get('dual_obstruction');parent_receipt=replay(parent)
    oracle=SpinZeroOracle(parent);states=parent['states']
    a=[[oracle.action(s).get(t,F(0)) for t in states] for s in states]
    result=separator_pairing(a,parent,certificate.get('indices'),certificate.get('signs'))
    result.update(states=[states[i] for i in certificate['indices']],
        parent_obstruction_source_states=parent_receipt['unique_action_states'],
        scope='The signed flat square in WHW coordinates pulls back to a positive rank-one atom in the physical coordinates and has negative exact pairing with the parent restricted dual. The normalized pairing uses unit physical trace for both atom and dual. This rejects that dual after adding the atom; it does not prove the enlarged cone sufficient. Complete replay also verifies the parent obstruction.')
    return result


def replay(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    if certificate.get('kind')!='spin_signed_clique_dual_v1':
        raise ValueError('Unsupported signed-edge clique dual')
    oracle=SpinZeroOracle(certificate);p=oracle.retained(certificate.get('retained_states'))
    states=certificate.get('states')
    if (type(states) is not list or not 1<=len(states)<=256
            or any(not oracle.valid_state(s) or s in p for s in states)
            or states!=sorted(set(states))):
        raise ValueError('Bounded distinct ordered physical Q support required')
    if certificate.get('metric_weights') is None: raise ValueError('Explicit diagonal metric required')
    a=[[oracle.action(s).get(t,F(0)) for t in states] for s in states]
    receipt=verify_signed(a,certificate)
    upper=F(receipt['packing_threshold_upper']);target=F(certificate['target_lower'])
    receipt.update(target_lower=rational_text(target),target_excluded=upper<target,
        unique_action_states=len(oracle.cache),referenced_determinants=oracle.referenced_state_count(),
        scope='Exact dual upper bound for balanced support-three/four clique squares and arbitrary real diagonally dominant residuals in the declared metric. Both signs of every edge square are allowed, so clique coefficients may overspend Hamiltonian edges and cancel. The dual is zero outside its declared Q support; restrictions of global balanced cliques to that support are checked by the lower-support inequalities. No complete sector or component coverage is needed. Does not exclude other signed three/four-coordinate atoms, arbitrary local amplitudes, or changed local metrics.')
    return receipt


def repair(a,weights,yn,zn,scale=10**9):
    """Round, clip pair inequalities, then add an exact multiple of the identity."""
    transformed,m,edges,b,groups=geometry(a,weights)
    if len(yn)!=len(a) or len(zn)!=len(edges): raise ValueError('Proposal dimensions do not match')
    y=[max(0,round(float(v)*scale)) for v in yn]
    if not any(y): y=[1]*len(a)
    z={e:min(2*(y[e[0]]+y[e[1]]),max(0,round(float(v)*scale))) for e,v in zip(edges,zn)}
    shift=0
    for group in groups:
        deficit=(len(group)-2)*sum(y[i] for i in group)-sum(z[e] for e in combinations(group,2))
        shift=max(shift,(deficit+len(group)-1)//len(group))
    # y_i+=t,z_ij+=2t adds tI to B. Clique slack grows by |C|t;
    # opposite-edge slack grows by 2t. Previously valid gates stay valid.
    y=[v+shift for v in y];z=[v+2*shift for v in z.values()]
    item={'metric_weights':weights,'node_weights':y,'edge_weights':z}
    return item,dict(verify_signed(a,item),identity_repair_integer=shift,rounding_scale=scale)


def construct(source,probe,output):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    source,probe,out=Path(source),Path(probe),Path(output)
    if out.exists(): raise ValueError('Preserve previous signed-clique export')
    old=json.loads(source.read_text());proposal=json.loads(probe.read_text())
    if old.get('kind')!='spin_clique_metric_probe_v1': raise ValueError('Metric probe required')
    index=proposal.get('block_index')
    if type(index) is not int or not 0<=index<len(old['blocks']): raise ValueError('Valid proposed block required')
    block=old['blocks'][index];states=block['states'];oracle=SpinZeroOracle(old)
    a=[[oracle.action(s).get(t,F(0)) for t in states] for s in states]
    item,diagnostic=repair(a,block['metric_weights'],proposal['node_weights_numerical'],proposal['edge_weights_numerical'])
    c={k:old[k] for k in ('modes','particles','hamiltonian','retained_states','target_lower')}
    c.update(kind='spin_signed_clique_dual_v1',states=states,**item)
    receipt=replay(c);out.mkdir(parents=True)
    diagnostic.update(source=str(source),numerical_proposal=str(probe))
    for name,data in [('certificate.json',c),('receipt.json',receipt),('proposal.json',diagnostic)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify');p.add_argument('--source')
    p.add_argument('--verify-separator')
    p.add_argument('--probe');p.add_argument('--output');args=p.parse_args()
    if args.verify: print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
    elif args.verify_separator: print(json.dumps(replay_separator(json.loads(Path(args.verify_separator).read_text())),indent=2))
    elif args.source and args.probe and args.output: print(json.dumps(construct(args.source,args.probe,args.output),indent=2))
    else: p.error('Specify --verify or --source, --probe and --output')
