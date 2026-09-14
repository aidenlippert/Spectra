"""Exact overlapping anticommuting certificates with untrusted float proposals.

Each checked split reconstructs the requested operator exactly. Optimization
convergence is unnecessary for soundness. No learned algorithm is claimed.
"""
from fractions import Fraction as F
from itertools import combinations
from math import sqrt, isfinite
from time import perf_counter
from .v7_certificate import proof_fraction, clean, norm_witness
from .certificates import _anti, _sqrt_interval

MAX_TERMS = 512
MAX_INCIDENCES = 2048
MAX_GROUPS = 512  # Ordinary singleton partitions must retain the V7 term budget.
MAX_PROPOSAL_GROUPS = 128


def _terms(op):
    if not isinstance(op, dict) or len(op) > MAX_TERMS:
        raise ValueError('term budget')
    if not op:
        return {}
    p = next(iter(op))
    if not isinstance(p, str) or not 1 <= len(p) <= 32:
        raise ValueError('Pauli width')
    return clean({p: proof_fraction(c) for p, c in op.items()}, len(p), MAX_TERMS)


def check_fractional_cover(op, witness):
    """Independent exact reconstruction, graph, and root-inequality checks."""
    start = perf_counter()
    cost = dict(pair_checks=0, coefficient_adds=0, squares=0)
    try:
        op = _terms(op)
        if not isinstance(witness, dict) or witness.get('schema') != 'v8-cover-2':
            raise ValueError('schema')
        groups = witness['groups']
        if not isinstance(groups, list) or len(groups) > MAX_GROUPS:
            raise ValueError('group budget')
        incidence = 0
        reconstruction = {p: F(0) for p in op}
        total = F(0)
        for group in groups:
            if not isinstance(group, dict) or set(group) != {'coefficients', 'upper'}:
                raise ValueError('group schema')
            coefficients = group['coefficients']
            if not isinstance(coefficients, dict) or not coefficients:
                raise ValueError('empty or invalid coefficients')
            incidence += len(coefficients)
            if incidence > MAX_INCIDENCES:
                raise ValueError('incidence budget')
            labels = list(coefficients)
            if any(p not in op for p in labels):
                raise ValueError('unknown term')
            for p, q in combinations(labels, 2):
                cost['pair_checks'] += 1
                if not _anti(p, q):
                    raise ValueError('non-anticommuting group')
            sq = F(0)
            for p, raw in coefficients.items():
                c = proof_fraction(raw)
                if not c:
                    raise ValueError('zero incidence')
                reconstruction[p] += c
                sq += c*c
                cost['coefficient_adds'] += 1
                cost['squares'] += 1
            upper = proof_fraction(group['upper'])
            cost['squares'] += 1
            if upper < 0 or upper*upper < sq:
                raise ValueError('invalid root upper bound')
            total += upper
        if reconstruction != op:
            raise ValueError('reconstruction mismatch')
        if proof_fraction(witness['claimed_bound']) != total:
            raise ValueError('claimed bound mismatch')
        return dict(status='certified', bound=str(total), checking_cost=cost,
                    checking_seconds=perf_counter()-start)
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, ZeroDivisionError, OverflowError) as exc:
        return dict(status='rejected', reason=str(exc), checking_cost=cost,
                    checking_seconds=perf_counter()-start)


def _witness(splits):
    groups = []
    total = F(0)
    for split in splits:
        split = {p: c for p, c in split.items() if c}
        if not split:
            continue
        square = sum((c*c for c in split.values()), F(0))
        upper = _sqrt_interval(square, 16)[1]
        total += upper
        groups.append(dict(coefficients={p: str(c) for p,c in split.items()}, upper=str(upper)))
    return dict(schema='v8-cover-2', groups=groups, claimed_bound=str(total))


def _group_family(op, supplied=None):
    """Union ordinary partitions and bounded seeded greedy *cliques*."""
    labels = sorted(op)
    groups = []
    known = set()
    incidence = 0
    pair_checks = 0
    def retain(group):
        nonlocal incidence
        group = tuple(sorted(group))
        if group in known:
            return
        if len(groups) >= MAX_PROPOSAL_GROUPS or incidence + len(group) > MAX_INCIDENCES:
            return
        groups.append(group)
        known.add(group)
        incidence += len(group)
    if supplied is not None:
        if not isinstance(supplied, (tuple,list)) or len(supplied)>MAX_PROPOSAL_GROUPS:
            raise ValueError('group budget')
        if sum(len(g) for g in supplied)>MAX_INCIDENCES:
            raise ValueError('incidence budget')
        for g in supplied:
            if not isinstance(g,(tuple,list)) or not g or len(set(g))!=len(g):
                raise ValueError('invalid group')
            if any(p not in op for p in g):
                raise ValueError('unknown term')
            for p,q in combinations(g,2):
                pair_checks += 1
                if not _anti(p,q):
                    raise ValueError('non-anticommuting group')
            retain(g)
    else:
        for mode in ('weighted','firstfit'):
            partition, stats = norm_witness(op,mode)
            pair_checks += stats['group_comparisons']
            for g in partition:
                retain(g['labels'])
        orders = [labels, sorted(op,key=lambda p:(-abs(op[p]),p)), list(reversed(labels))]
        for order in orders:
            for seed in order:
                if len(groups)>=MAX_PROPOSAL_GROUPS or incidence>=MAX_INCIDENCES:
                    break
                group=[seed]
                for p in order:
                    if p == seed:
                        continue
                    compatible=True
                    for q in group:
                        pair_checks+=1
                        if not _anti(p,q):
                            compatible=False
                            break
                    if compatible:
                        group.append(p)
                retain(group)
    if set(p for g in groups for p in g) != set(op):
        raise ValueError('uncovered terms within group budget')
    return groups, dict(group_pair_checks=pair_checks, group_count=len(groups), incidences=incidence)


def propose_fractional_cover(op, *, groups=None, iterations=16):
    """Bounded cyclic coordinate minimization; return best checked feasible split.

For one coefficient c, minimizing sum_j sqrt(a_j^2+x_j^2), sum_j x_j=c,
sets x_j=c*a_j/sum(a). Here a_j is the norm of the other entries in group j.
The formula follows from the triangle inequality on vectors (a_j,x_j).
All floats are untrusted; exact coefficient reconciliation follows optimization.
"""
    start=perf_counter()
    op=_terms(op)
    if type(iterations) is not int or not 0<=iterations<=128:
        raise ValueError('iteration budget')
    if not op:
        w=_witness([])
        return w,dict(construction_seconds=perf_counter()-start, iterations=0)
    family,stats=_group_family(op,groups)
    inc={p:[j for j,g in enumerate(family) if p in g] for p in op}
    # Scale to avoid overflow in the untrusted optimizer; exact bounds have no scale restriction.
    scale=max(abs(c) for c in op.values())
    vals={p:float(c/scale) for p,c in op.items()}
    x=[{p:vals[p]/len(inc[p]) for p in g} for g in family]
    def rationalized():
        splits=[{p: F.from_float(v).limit_denominator(10**9)*scale for p,v in row.items()} for row in x]
        for p in op:
            residue=op[p]-sum((splits[j][p] for j in inc[p]), F(0))
            # Reconcile into largest incidence to avoid an avoidable large relative correction.
            j=max(inc[p],key=lambda j:abs(splits[j][p]))
            splits[j][p]+=residue
        return _witness(splits)
    best=rationalized()
    # Conventional partitions remain eligible: extra search must never force a looser bound.
    for mode in ('l1','firstfit','weighted'):
        partition,cost=norm_witness(op,mode)
        stats['group_pair_checks']+=cost['group_comparisons']
        if len(partition)<=MAX_GROUPS:
            w=_witness([{p:op[p] for p in g['labels']} for g in partition])
            if F(w['claimed_bound'])<F(best['claimed_bound']):
                best=w
    updates=0
    objective=[]
    for sweep in range(iterations):
        for p in sorted(op):
            js=inc[p]
            if len(js)<2:
                continue
            a=[sqrt(sum(v*v for q,v in x[j].items() if q!=p)) for j in js]
            total=sum(a)
            if total:
                for j,v in zip(js,a):
                    x[j][p]=vals[p]*v/total
                    updates+=1
        value=sum(sqrt(sum(v*v for v in row.values())) for row in x)
        if not isfinite(value):
            raise ValueError('nonfinite proposal')
        objective.append(value)
    w=rationalized()
    if F(w['claimed_bound'])<F(best['claimed_bound']):
        best=w
    receipt=check_fractional_cover(op,best)
    if receipt['status']!='certified':
        raise ValueError('proposed cover failed exact verification: '+receipt['reason'])
    stats.update(iterations=iterations, coordinate_updates=updates, scaled_objective_trace=objective,
                 construction_seconds=perf_counter()-start, internal_check=receipt)
    return best,stats


def check_evolution_cover(gen, initial, pieces, witness, tolerance, *, expected_time):
    """Same exact residual recomputation for partition and overlapping-cover arms."""
    from .v7_certificate import rational, residual_records
    start=perf_counter()
    try:
        tolerance=rational(tolerance); expected_time=rational(expected_time)
        if tolerance<0 or expected_time<=0:
            raise ValueError('invalid request')
        if witness['schema']!='v8-evolution-cover-1':
            raise ValueError('schema')
        records=residual_records(gen,initial,pieces,'power')
        if sum((p.duration for p in pieces),F(0))!=expected_time:
            raise ValueError('horizon mismatch')
        if set(witness['witnesses'])!={label for label,_,_ in records}:
            raise ValueError('residual coverage')
        total=F(0); checks=[]
        for label,op,weight in records:
            r=check_fractional_cover(op,witness['witnesses'][label])
            checks.append(r)
            if r['status']!='certified':
                raise ValueError(label+': '+r['reason'])
            total+=weight*F(r['bound'])
        if total!=proof_fraction(witness['claimed_bound']):
            raise ValueError('claimed mismatch')
        return dict(status='certified' if total<=tolerance else 'over_tolerance',bound=str(total),
                    norm_checks=checks,generator_cost=dict(gen.cost),checking_seconds=perf_counter()-start)
    except (ValueError,TypeError,KeyError,IndexError,AttributeError,ZeroDivisionError,OverflowError) as exc:
        return dict(status='rejected',reason=str(exc),generator_cost=dict(gen.cost),checking_seconds=perf_counter()-start)


def adaptive_cover(gen,initial,duration,tolerance,*,enabled=False,max_order=24):
    """Fixed prospective trigger: failed ordinary bound <=1.1 eps, <=40 terms."""
    from .v7_certificate import Piece, rational
    duration=rational(duration); tolerance=rational(tolerance)
    if duration<=0 or tolerance<0 or type(max_order) is not int or not 0<=max_order<=24:
        raise ValueError('request/order')
    coeffs=[clean(initial,gen.n,gen.max_terms)]; probes=[]; cover_calls=[]
    empty=_witness([])
    for m in range(max_order+1):
        derivative=gen.apply(coeffs[-1]);factor=duration**(m+1)/F(m+1)
        best=None;bound=None
        for mode in ('l1','firstfit','weighted'):
            groups,_=norm_witness(derivative,mode)
            w=dict(schema='v8-cover-2', groups=[dict(coefficients={p:str(derivative[p]) for p in g['labels']},upper=g['upper']) for g in groups],
                   claimed_bound=str(sum((F(g['upper']) for g in groups),F(0))))
            b=F(w['claimed_bound'])*factor
            probes.append(dict(order=m,mode=mode,bound=str(b)))
            if bound is None or b<bound:bound=b;best=w
            if bound<=tolerance:break
        if enabled and tolerance<bound<=F(11,10)*tolerance and len(derivative)<=40:
            candidate,cost=propose_fractional_cover(derivative,iterations=16)
            b=F(candidate['claimed_bound'])*factor
            cover_calls.append(dict(order=m,bound=str(b),cost=cost))
            if b<bound:bound=b;best=candidate
        if bound<=tolerance:
            ws={'jump:0':empty,**{f'residual:0:{k}':empty for k in range(m+1)}}
            if derivative:
                # Residual is -derivative: negation preserves all group norms.
                best=dict(best,groups=[dict(coefficients={p:str(-F(c)) for p,c in g['coefficients'].items()},upper=g['upper']) for g in best['groups']])
                ws[f'residual:0:{m}']=best
            else:ws={'jump:0':empty,'residual:0:0':empty}
            return Piece(duration,tuple(coeffs)),dict(schema='v8-evolution-cover-1',witnesses=ws,claimed_bound=str(bound)),dict(probes=probes,cover_calls=cover_calls,generator_cost=dict(gen.cost))
        if m<max_order:coeffs.append({p:c/F(m+1) for p,c in derivative.items()})
    raise ValueError('Taylor order budget')
