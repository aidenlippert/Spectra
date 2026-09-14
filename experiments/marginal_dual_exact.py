"""Exact feasible cubic moment witnesses, with optional numerical discovery."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
import argparse
import json
from pathlib import Path

from experiments.marginal_symbolic import adj, canonical, decode, encode, mono, multiplier_basis, number_shift, product
from experiments.marginal_collective import hopping_polynomial


def component_charge(w,modes):
    m=modes//2;result=[0]*m
    for c,i in w:result[i%m]+=2*c-1
    return tuple(result)


def rows_for(modes):
    rows=[]
    for k in range(4):
        sets=list(combinations(range(modes),k))
        for i,left in enumerate(sets):
            for right in sets[i:]:
                w=tuple((1,j) for j in left)+tuple((0,j) for j in right)
                if not any(component_charge(w,modes)):rows.append(w)
    return rows


def full_word_blocks(modes,max_degree=3):
    groups={}
    for degree in range(max_degree+1):
        for nc in range(degree+1):
            for left in combinations(range(modes),nc):
                for right in combinations(range(modes),degree-nc):
                    w=tuple((1,i) for i in left)+tuple((0,i) for i in right)
                    groups.setdefault(component_charge(w,modes),[]).append(w)
    return [{'name':str(q),'words':words} for q,words in sorted(groups.items())]


def full_cubic_blocks(modes):
    return full_word_blocks(modes,3)


def functional(rows,values,modes):
    table={w:F(v) for w,v in zip(rows,values)}
    def evaluate(poly):
        answer=F(0)
        for w,c in canonical(poly).items():
            if any(component_charge(w,modes)):continue
            left=tuple(i for cr,i in w if cr);right=tuple(i for cr,i in w if not cr)
            if len(left)!=len(right):raise ValueError('Nonconserving functional query')
            key=w if left<=right else tuple((1,i) for i in right)+tuple((0,i) for i in left)
            answer+=c*table[key]/(1 if left==right else 2)
        return answer
    return evaluate


def exact_psd(matrix):
    """Exact symmetric LDL elimination, including mandatory zero pivots."""
    n=len(matrix);a=[list(map(F,row)) for row in matrix];rank=0
    if any(len(row)!=n for row in a) or any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Matrix must be symmetric')
    for k in range(n):
        pivot=a[k][k]
        if pivot<0:return False,rank
        if pivot==0:
            if any(a[j][k] for j in range(k+1,n)):return False,rank
            continue
        rank+=1
        for i in range(k+1,n):
            if not a[i][k]:continue
            ratio=a[i][k]/pivot
            for j in range(i,n):
                a[j][i]-=ratio*a[j][k]
                a[i][j]=a[j][i]
    return True,rank


def verify(witness):
    modes=witness['modes'];n=witness['particles']
    if type(modes) is not int or modes<4 or modes%2 or type(n) is not int or n!=modes//2:
        raise ValueError('Expected matched half-filled sector')
    if not isinstance(witness['t'],str):raise ValueError('Exact t string required')
    h=decode(witness['hamiltonian'],modes,4)
    if h!=hopping_polynomial(modes,F(witness['t'])):raise ValueError('Hamiltonian mismatch')
    rows=rows_for(modes);values=witness['moments']
    if len(values)!=len(rows) or any(not isinstance(v,str) for v in values):raise ValueError('Invalid rational moments')
    evaluate=functional(rows,values,modes)
    if evaluate(mono(()))!=1:raise ValueError('Moment normalization failed')
    shift=number_shift(modes,n);equalities=0
    for p in multiplier_basis(modes,2):
        if evaluate(product(shift,p))!=0:raise ValueError('Number equality failed')
        equalities+=1
    ranks=[]
    for block in full_cubic_blocks(modes):
        words=block['words']
        matrix=[[evaluate(product(adj(mono(u)),mono(v))) for v in words] for u in words]
        valid,rank=exact_psd(matrix)
        if not valid:raise ValueError(f"Non-PSD moment block {block['name']}")
        ranks.append(rank)
    energy=evaluate(h)
    return {'energy':str(energy),'energy_float':float(energy),'moment_coordinates':len(rows),
            'number_equalities_checked':equalities,'psd_blocks':len(ranks),'block_ranks':ranks}


def build(modes=10,t=F(1,5)):
    import numpy as np
    import sympy as sy
    from experiments.marginal_adaptive import assemble_and_solve
    h=hopping_polynomial(modes,t);blocks=full_cubic_blocks(modes)
    proposal,data=assemble_and_solve(h,modes,modes//2,blocks,groups=[[i,i+modes//2] for i in range(modes//2)])
    rows=data['rows']; assert rows==rows_for(modes)
    y=data['dual'];a=data['multiplier_map'].T.toarray()
    # All map coefficients are integer or half-integer, reconstructed exactly.
    matrix=sy.Matrix([[sy.Rational(int(round(2*v)),2) for v in row]+[0] for row in a])
    unit=[0]*len(rows);unit[rows.index(())]=1
    matrix=matrix.col_join(sy.Matrix([unit+[1]]))
    reduced,pivots=matrix.rref();assert len(rows) not in pivots
    free=[j for j in range(len(rows)) if j not in pivots]
    values=[F(0)]*len(rows)
    for j in free:values[j]=F(int(round(y[j]*10**10)),10**10)
    for i,j in enumerate(pivots):
        values[j]=F(reduced[i,-1])-sum(F(reduced[i,k])*values[k] for k in free)
    physical=[]
    for w in rows:
        k=len(w)//2;left=tuple(i for cr,i in w if cr);right=tuple(i for cr,i in w if not cr)
        physical.append(F((-1)**(k*(k-1)//2)*comb(modes//2,k),comb(modes,k)) if left==right else F(0))
    out=Path(__file__).resolve().parents[1]/'results/marginal_dual';out.mkdir(exist_ok=True)
    (out/'numeric.json').write_text(json.dumps({'energy':float(data['rhs']@y),'normalization':float(y[rows.index(())]),
         'equality_residual':float(np.max(np.abs(a@y))),'status':proposal['status'],'affine_rank':len(pivots),'coordinates':len(rows)},indent=2)+'\n')
    attempts=[]
    for mixing in (F(1,100000),F(1,10000),F(1,1000),F(1,100)):
        repaired=[(1-mixing)*v+mixing*p for v,p in zip(values,physical)]
        witness={'modes':modes,'particles':modes//2,'t':str(t),'hamiltonian':encode(h),'moments':list(map(str,repaired)),
                 'physical_mixing':str(mixing)}
        try:receipt=verify(witness)
        except ValueError as error:
            attempts.append({'mixing':str(mixing),'error':str(error)});continue
        (out/'witness.json').write_text(json.dumps(witness)+'\n')
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps(receipt),flush=True)
        (out/'attempts.json').write_text(json.dumps(attempts,indent=2)+'\n')
        return witness,receipt
    (out/'attempts.json').write_text(json.dumps(attempts,indent=2)+'\n')
    raise RuntimeError('No exact PSD repair; see attempts.json')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--build',action='store_true');parser.add_argument('--verify',type=Path)
    args=parser.parse_args()
    if args.build:build()
    elif args.verify:print(json.dumps(verify(json.loads(args.verify.read_text())),indent=2))
    else:parser.error('Choose --build or --verify')
