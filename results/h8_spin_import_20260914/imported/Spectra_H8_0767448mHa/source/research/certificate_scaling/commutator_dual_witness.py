"""Exact dual obstruction for an explicitly rebuilt H-only commutator cone.

The verifier imports only the standard library and exact CAR algebra. It
rebuilds the full charge blocks (without symmetry pruning), every body-two
number-ideal constraint, the coefficient box and all rational PSD checks.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import argparse, hashlib, json, sys, time
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import (canonical, hermitian, decode, encode,
    add, scale, mono, product, number_shift, multiplier_basis, verify as primal_verify)
from experiments.marginal_hunt_car import adj


def groups(h, m, enrich, channels):
    minus = [mono(((0, i),)) for i in range(m)]
    pairs = [canonical(mono(((0, j), (0, i)))) for i, j in combinations(range(m), 2)]
    number = [mono(((1, i), (0, j))) for i in range(m) for j in range(m)]
    if enrich:
        for i in range(m):
            a = mono(((0, i),))
            cubic = {w:c for w,c in add(product(h,a), scale(product(a,h),-1)).items() if len(w)==3}
            if channels:
                minus.extend({w:c for w,c in cubic.items() if w[0][1]==j}
                             for j in sorted({w[0][1] for w in cubic}))
            elif cubic:
                minus.append(cubic)
    return [minus, [canonical(adj(p)) for p in minus], pairs,
            [canonical(adj(p)) for p in pairs], number]


def evaluate(p, moments):
    return sum((c*moments.get(w,F(0)) for w,c in p.items()), F(0))


def gram(polys, moments):
    # Full matrix entries y(p_i^dagger p_j). Off-diagonal symmetrized
    # upper-triangle primal columns would have to be divided by TWO here.
    matrix = [[F(0) for _ in polys] for _ in polys]
    for i,p in enumerate(polys):
        for j in range(i,len(polys)):
            value=evaluate(product(canonical(adj(p)),polys[j]),moments)
            reverse=evaluate(product(canonical(adj(polys[j])),p),moments)
            if value!=reverse:raise ValueError('Non-Hermitian moment Gram')
            matrix[i][j]=matrix[j][i]=value
    return matrix


def psd(matrix):
    """Exact symmetric elimination, with necessary zero-pivot row refusal."""
    if any(type(v) is not int and not isinstance(v,F) for row in matrix for v in row):
        raise ValueError('Exact rational matrix entries required')
    a=[[F(v) for v in row] for row in matrix]; n=len(a);rank=0
    if any(len(row)!=n for row in a) or any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Non-symmetric square matrix')
    for k in range(n):
        pivot=a[k][k]
        if pivot<0:raise ValueError('Negative exact PSD pivot')
        if not pivot:
            if any(a[k][j] for j in range(k+1,n)):raise ValueError('Nonzero row at zero PSD pivot')
            continue
        rank+=1
        for i in range(k+1,n):
            if not a[i][k]:continue
            ratio=a[i][k]/pivot
            for j in range(i,n):
                a[i][j]-=ratio*a[k][j];a[j][i]=a[i][j]
    return {'dimension':n,'rank':rank}


def moment_decode(items,m):
    # Moments are a functional, not a polynomial: never canonicalize them
    # by adding coefficients as though they represented an operator.
    result={}
    from experiments.marginal_symbolic import validate_word
    for item in items:
        w=validate_word(item['word'],m,6)
        if canonical(mono(w))!=mono(w) or w in result:raise ValueError('Duplicate/noncanonical moment')
        if not isinstance(item['value'],str):raise ValueError('Exact rational moment required')
        result[w]=F(item['value'])
    return result


def check(witness, comparison=None):
    start=time.monotonic();m=witness['modes'];n=witness['particles']
    if witness.get('one_body_steps',0)!=0:raise ValueError('Krylov-enriched dictionaries are outside this checker')
    if type(m) is not int or m<4 or type(n) is not int or not 0<n<m:raise ValueError('Invalid sector')
    if type(witness['enrich']) is not bool or type(witness['creator_channels']) is not bool:raise ValueError('Invalid dictionary flags')
    h=decode(witness['hamiltonian'],m,4)
    if not hermitian(h) or any(sum(2*c-1 for c,_ in w) for w in h):raise ValueError('Invalid Hamiltonian')
    y=moment_decode(witness['moments'],m)
    if y.get(())!=1 or any(abs(v)>1 for v in y.values()):raise ValueError('Normalization/coefficient box failed')
    for w,v in y.items():
        if evaluate(canonical(adj(mono(w))),y)!=v:raise ValueError('Hermitian functional failed')
    ideal=[product(number_shift(m,n),p) for p in multiplier_basis(m,max_body=2)]
    if any(evaluate(p,y) for p in ideal):raise ValueError('Number-ideal dual equality failed')
    stats=[psd(gram(g,y)) for g in groups(h,m,witness['enrich'],witness['creator_channels'])]
    objective=evaluate(h,y)
    receipt={'dual_objective':str(objective),'dual_objective_float':float(objective),
             'full_charge_Gram_checks':stats,'ideal_equalities':len(ideal),
             'moment_count':len(y),'max_moment_bits':max(max(v.numerator.bit_length(),v.denominator.bit_length()) for v in y.values()),
             'scope':'Upper bound on every b-||residual||_1 from this rebuilt commutator dictionary and body-two number ideal; not a ground-energy lower bound.'}
    if comparison is not None:
        if comparison['modes']!=m or comparison['particles']!=n or decode(comparison['hamiltonian'],m,4)!=h:
            raise ValueError('Comparison Hamiltonian/sector differs')
        physical=primal_verify(comparison);gap=F(physical['lower'])-objective
        receipt.update({'physical_lower':physical['lower'],'proved_accuracy_obstruction':str(gap),
                        'proved_accuracy_obstruction_float':float(gap),'excludes_0_0016':gap>F(16,10000)})
    receipt['replay_seconds']=time.monotonic()-start
    return receipt


def seed(w,m,n):
    if not w:return F(1)
    left=tuple(i for c,i in w if c);right=tuple(i for c,i in w if not c);k=len(left)
    return F((-1)**(k*(k-1)//2)*comb(n,k),comb(m,k)) if left==right and k<=n else F(0)


def propose(raw, out, comparison):
    """Rational affine projection, then exact-tested mixing with sector trace."""
    if raw.get('one_body_steps',0)!=0:raise ValueError('Krylov-enriched dictionaries are outside this checker')
    start=time.monotonic();m=raw['modes'];n=raw['particles'];h=decode(raw['hamiltonian'],m,4)
    proposed={tuple(tuple(x) for x in w):F(round(v*10**12),10**12) for w,v in zip(raw['rows'],raw['dual'])}
    allwords=set(proposed)
    ideals=[mono(())]+[product(number_shift(m,n),p) for p in multiplier_basis(m,max_body=2)]
    for p in ideals:allwords.update(p)
    # Only store moments in the numerical symmetry support. Missing moments
    # are zero; all full (unpruned) ideal equalities are nevertheless checked.
    active=set(proposed)
    representatives={}
    for w in active:
        aw=canonical(adj(mono(w)))
        if len(aw)!=1:raise ValueError('Unexpected adjoint expansion')
        partner,sign=next(iter(aw.items()))
        if sign!=1:raise ValueError('Balanced normal-order adjoint sign changed')
        representatives[w]=min(w,partner)
    keys=sorted(set(representatives.values()),key=lambda w:(len(w),w));lookup={w:i for i,w in enumerate(keys)}
    values=[]
    for w in keys:
        partner=next(iter(canonical(adj(mono(w)))))
        values.append((proposed.get(w,F(0))+proposed.get(partner,F(0)))/2)
    # Sparse exact forward elimination on integral ideal constraints.
    pivots={}
    for rownum,p in enumerate(ideals):
        row={}
        for w,c in p.items():
            if w in representatives:
                j=lookup[representatives[w]];row[j]=row.get(j,F(0))+c
        row={j:v for j,v in row.items() if v};rhs=F(rownum==0)
        while row:
            j=min(row)
            if j not in pivots:
                v=row[j];pivots[j]=({k:c/v for k,c in row.items()},rhs/v);break
            old,b=pivots[j];v=row[j];rhs-=v*b
            for k,c in old.items():
                row[k]=row.get(k,F(0))-v*c
                if not row[k]:del row[k]
        else:
            if rhs:raise ValueError('Inconsistent dual affine constraints')
    for j,(row,rhs) in sorted(pivots.items(),reverse=True):
        values[j]=rhs-sum(c*values[k] for k,c in row.items() if k!=j)
    rounded={w:values[lookup[r]] for w,r in representatives.items()}
    # Fill all balanced degree<=6 diagonal moments for the trace seed. The
    # affine candidate has zero outside its reported support.
    for k in range(4):
        for inds in combinations(range(m),k):allwords.add(tuple((1,i) for i in inds)+tuple((0,i) for i in inds))
    attempts=[];out=Path(out);out.mkdir(parents=True,exist_ok=False)
    for mix in (F(0),F(1,10**8),F(1,10**7),F(1,10**6),F(1,10**5),F(1,10**4),F(1,1000)):
        y={w:(1-mix)*rounded.get(w,F(0))+mix*seed(w,m,n) for w in allwords}
        witness={'modes':m,'particles':n,'hamiltonian':encode(h),'enrich':raw['enrich'],
                 'creator_channels':raw['creator_channels'],'trace_mixture':str(mix),
                 'moments':[{'word':[list(x) for x in w],'value':str(v)} for w,v in sorted(y.items(),key=lambda x:(len(x[0]),x[0])) if v or not w]}
        try:receipt=check(witness,comparison)
        except ValueError as e:
            attempts.append({'mix':str(mix),'refusal':str(e)});continue
        receipt.update({'attempts':attempts,'trace_mixture':str(mix),'affine_rank':len(pivots),'construction_seconds':time.monotonic()-start})
        text=json.dumps(witness,separators=(',',':'))+'\n';(out/'witness.json').write_text(text)
        receipt['witness_sha256']=hashlib.sha256(text.encode()).hexdigest();receipt['witness_bytes']=len(text.encode())
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt
    (out/'refusals.json').write_text(json.dumps(attempts,indent=2)+'\n')
    raise ValueError('No exact positive dual repair accepted')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--witness',type=Path);p.add_argument('--proposal',type=Path)
    p.add_argument('--comparison',type=Path,required=True);p.add_argument('--out',type=Path)
    a=p.parse_args();comparison=json.loads(a.comparison.read_text())
    if bool(a.witness)==bool(a.proposal):p.error('Exactly one witness or proposal required')
    result=(check(json.loads(a.witness.read_text()),comparison) if a.witness else propose(json.loads(a.proposal.read_text()),a.out,comparison))
    print(json.dumps(result,indent=2))
