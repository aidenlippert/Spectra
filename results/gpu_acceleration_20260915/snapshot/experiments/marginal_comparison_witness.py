"""Hamiltonian-bound exact FW-k comparison duals; never physical lower bounds."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
from experiments.marginal_spin_reduction import SpinZeroOracle


def _check(h, v, k):
    n=len(h)
    if type(k) is not int or k not in (2,3,4): raise ValueError('k must be 2, 3, or 4')
    if not 1<=n<=512 or len(v)!=n or any(type(x) is not int or not 0<x<=10**12 for x in v):
        raise ValueError('Bounded nonempty positive integer vector required')
    if any(len(row)!=n for row in h) or any(h[i][j]!=h[j][i] for i in range(n) for j in range(i)):
        raise ValueError('Symmetric square matrix required')


def witness_matrix(h, v, k):
    _check(h,v,k)
    return [[F(v[i]*v[i]) if i==j else F(-((h[i][j]>0)-(h[i][j]<0))*v[i]*v[j],k-1)
             for j in range(len(v))] for i in range(len(v))]


def ceiling(h, v, k):
    _check(h,v,k)
    diagonal=sum((F(h[i][i])*v[i]*v[i] for i in range(len(v))),F(0))
    edges=sum((abs(F(h[i][j]))*v[i]*v[j] for i,j in combinations(range(len(v)),2)),F(0))
    return (diagonal-F(2,k-1)*edges)/sum(x*x for x in v)


def principal_psd(x, k):
    """Small test oracle only; production proof uses diagonal congruence."""
    if len(x)>8 or type(k) is not int or not 1<=k<=4: raise ValueError('Small principal-minor test required')
    def det(a):
        if len(a)==1: return a[0][0]
        return sum(((-1)**j*a[0][j]*det([row[:j]+row[j+1:] for row in a[1:]]) for j in range(len(a))),F(0))
    return all(det([[F(x[i][j]) for j in ids] for i in ids])>=0
               for size in range(1,k+1) for ids in combinations(range(len(x)),size))


def physical_matrix(oracle, retained_states, states):
    retained=oracle.retained(retained_states)
    if (type(states) is not list or not 1<=len(states)<=512
            or any(not oracle.valid_state(s) or s in retained for s in states)
            or len(set(states))!=len(states)):
        raise ValueError('Distinct physical complement witness states required')
    columns=[oracle.action(s) for s in states]
    return [[columns[j].get(s,F(0)) for j in range(len(states))] for s in states]


def replay(data):
    if data.get('kind')!='spin_comparison_family_dual_v1': raise ValueError('Unsupported physical comparison dual')
    oracle=SpinZeroOracle(data)
    h=physical_matrix(oracle,data.get('retained_states'),data.get('witness_states'))
    v,k=data.get('v'),data.get('k'); got=ceiling(h,v,k)
    if F(data['ceiling'])!=got: raise ValueError('Declared ceiling does not match the physical Hamiltonian')
    target=F(data['target_lower'])
    # On any <=k support, D_v^-1 X D_v has diagonal 1 and
    # absolute off-diagonal row sum <=(k-1)/(k-1)=1. Padding
    # outside witness_states by zero preserves this PSD guarantee.
    return {'ceiling':str(got),'ceiling_float':float(got),'k':k,'witness_support':len(v),
            'target_lower':str(target),'target_excluded':got<target,
            'target_minus_ceiling':str(target-got),
            'unique_action_states':len(oracle.cache),'referenced_determinants':oracle.referenced_state_count(),
            'scope':'Exact trace-normalized witness PSD on every principal support of size at most k by diagonal congruence and diagonal dominance, embedded by zero in the physical determinant complement. Its trace pairing gives an upper ceiling on FW-k lower certificates for this H/P/basis. This is not a physical energy lower bound; failure to exclude the target does not prove feasibility.'}


def construct(source, retained_states, output, k, target):
    from experiments.marginal_spin_constructor import spin_blocks
    import numpy as np
    output=Path(output)
    if output.exists(): raise ValueError('Preserve previous physical comparison dual')
    data=json.loads(Path(source).read_text()); base={key:data[key] for key in ('modes','particles','hamiltonian')}
    oracle=SpinZeroOracle(base); oracle.retained(retained_states)
    best=None
    for group in spin_blocks(oracle,retained_states):
        h=physical_matrix(oracle,retained_states,group)
        a=np.array(h,dtype=float); c=-np.abs(a)/(k-1); np.fill_diagonal(c,np.diag(a))
        vec=np.linalg.eigh(c)[1][:,0]; v=[max(1,int(round(abs(x)*10**8))) for x in vec]
        value=ceiling(h,v,k)
        if best is None or value<best[0]: best=(value,group,v)
    if best is None: raise ValueError('Nonempty complement required')
    value,group,v=best
    cert=dict(base,kind='spin_comparison_family_dual_v1',retained_states=retained_states,
              witness_states=group,v=v,k=k,ceiling=str(value),target_lower=str(F(target)))
    receipt=replay(cert)
    output.mkdir(parents=True)
    (output/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--verify',required=True);args=p.parse_args()
    print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
