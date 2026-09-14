"""Exact joint hopping envelope for finite-range positive charge products.

For each fixed doublon count D, maximize the SUM of hopping penalties by
finite-memory dynamic programming. Individual transition maxima are never
combined. This verifier is specific to the neutral open uniform Hubbard chain.
The metric is g(D) times positive local factors, with g(0)=0 and g(D)>0 for
D>0. Omitting the count profile recovers g(D)=D.

For any charge pattern, alternate spins within every singly occupied run.
Even-length runs are balanced; odd-length runs occur in an even number because
the total number of single sites is even. Choose half their phases each way.
This realizes a balanced-spin assignment with every single-single hop allowed.
Other charge-edge hopping multiplicities are independent of these phases.
Thus the charge-row maximum is attainable, and also bounds every spin row.
"""
from fractions import Fraction as F
from itertools import product
import json

from experiments.marginal_symbolic import decode, mono, product as mul, scale, add


def model(certificate):
    if certificate.get('kind') != 'joint_charge_product_dp_v1':
        raise ValueError('Unsupported joint charge product certificate')
    sites = certificate.get('sites')
    if type(sites) is not int or not 2 <= sites <= 32 or sites % 2:
        raise ValueError('Even site count in2..32 required')
    if certificate.get('modes') != 2*sites or certificate.get('particles') != sites:
        raise ValueError('Neutral half-filled sector required')
    def rational(value, positive=False):
        if type(value) is not str or len(value) > 128:
            raise ValueError('Bounded exact rational string required')
        result = F(value)
        if max(abs(result.numerator), result.denominator) > 10**30 or positive and result <= 0:
            raise ValueError('Rational magnitude or positivity gate failed')
        return result
    U, t = rational(certificate['U']), rational(certificate['t'])
    gamma = rational(certificate['target_lower'])
    if U < 0 or t < 0:
        raise ValueError('Nonnegative U and t required')
    occupations = [mono(((1,i),(0,i))) for i in range(2*sites)]
    expected = add(*(scale(mul(occupations[2*i], occupations[2*i+1]), U) for i in range(sites)),
                   *(mono(((1,2*i+s),(0,2*(i+1)+s)), -t) for i in range(sites-1) for s in (0,1)),
                   *(mono(((1,2*(i+1)+s),(0,2*i+s)), -t) for i in range(sites-1) for s in (0,1)))
    if decode(certificate['hamiltonian'], 2*sites, 4) != expected:
        raise ValueError('Hamiltonian differs from declared uniform open Hubbard chain')
    onsite = certificate.get('onsite')
    if type(onsite) is not list or len(onsite) != sites:
        raise ValueError('One positive charge table per site required')
    def table(values, size):
        if type(values) is not list or len(values) != size:
            raise ValueError('Invalid factor table shape')
        return [rational(v, True) for v in values]
    onsite = [table(v,3) for v in onsite]
    raw = certificate.get('pairs')
    if type(raw) is not list or len(raw) > 3*sites:
        raise ValueError('Bounded pair factor list required')
    pairs, seen = [], set()
    for pair in raw:
        i,j = pair['i'],pair['j']
        if type(i) is not int or type(j) is not int or not 0 <= i < j < sites or j-i > 3 or (i,j) in seen:
            raise ValueError('Unique pair range in1..3 required')
        values = pair['values']
        if type(values) is not list or len(values) != 3:
            raise ValueError('Pair table must be3x3')
        pairs.append((i,j,[table(row,3) for row in values])); seen.add((i,j))
    profile=certificate.get('doublon_weights')
    if profile is None:profile=[F(i) for i in range(sites//2+1)]
    else:
        if type(profile) is not list or len(profile)!=sites//2+1:
            raise ValueError('One count weight per feasible doublon count required')
        if rational(profile[0])!=0:
            raise ValueError('Count metric must vanish on the valence space')
        profile=[F(0)]+[rational(value,True) for value in profile[1:]]
    return sites,U,t,gamma,onsite,pairs,profile


def joint_envelope(certificate, max_states=500000, max_local_entries=200000, return_witnesses=False):
    if type(return_witnesses) is not bool:
        raise ValueError('Boolean witness option required')
    if type(max_states) is not int or not 1 <= max_states <= 500000 or type(max_local_entries) is not int or not 1 <= max_local_entries <= 200000:
        raise ValueError('Bounded positive DP budgets required')
    sites,U,t,gamma,onsite,pairs,profile = model(certificate)
    radius = max((j-i for i,j,_ in pairs), default=0)
    memory = 2*radius+1
    edges, ending = [], [[] for _ in range(sites)]
    local_entries = 0;largest_window=0
    for edge in range(sites-1):
        touched = [(i,j,tab) for i,j,tab in pairs if i in (edge,edge+1) or j in (edge,edge+1)]
        support = {edge,edge+1}
        for i,j,_ in touched: support.update((i,j))
        lo,hi = min(support),max(support)
        largest_window=max(largest_window,hi-lo+1)
        local_entries += 3**(hi-lo+1)
        if local_entries > max_local_entries:
            raise ValueError('Local window table budget exceeded')
        # Separate count changes -1,0,+1; their profile coefficients depend
        # only on the fixed total D, not on the eliminated charge prefix.
        table = {}
        for values in product((-1,0,1), repeat=hi-lo+1):
            q = dict(zip(range(lo,hi+1),values)); a,b = q[edge],q[edge+1]
            targets=[]
            if abs(a-b)==1: targets=[(b,a,1,0)]
            elif {a,b}=={-1,1}: targets=[(0,0,2,-1)]
            elif a==b==0: targets=[(1,-1,1,1),(-1,1,1,1)]
            coefficients=[F(0),F(0),F(0)]
            for na,nb,rate,delta in targets:
                new=dict(q);new[edge],new[edge+1]=na,nb
                ratio=onsite[edge][na+1]/onsite[edge][a+1]*onsite[edge+1][nb+1]/onsite[edge+1][b+1]
                for i,j,tab in touched:
                    ratio *= tab[new[i]+1][new[j]+1]/tab[q[i]+1][q[j]+1]
                coefficients[delta+1] += rate*ratio
            table[values]=coefficients
        edges.append((lo,hi,table));ending[hi].append(edge)
    rows=[];peak=0;visited=0
    for D in range(1,sites//2+1):
        count_weights=[profile[D-1],profile[D],profile[D+1] if D+1<len(profile) else F(0)]
        # State: prefix charge, prefix squared charge, required recent charges.
        states={(0,0,()):F(0)}
        paths={(0,0,()):()} if return_witnesses else None
        for k in range(sites):
            next_states={};remaining=sites-k-1
            next_paths={} if return_witnesses else None
            for (charge,squares,tail),score in states.items():
                for q in (-1,0,1):
                    cq,sq=charge+q,squares+q*q
                    if abs(cq)>remaining or sq>2*D or sq+remaining<2*D:
                        continue
                    sequence=tail+(q,);base=k-len(tail);value=score
                    for edge in ending[k]:
                        lo,hi,table=edges[edge]
                        coefficients=table[sequence[lo-base:hi-base+1]]
                        value+=sum(a*b for a,b in zip(coefficients,count_weights))
                    key=(cq,sq,sequence[-memory:])
                    if key not in next_states or value>next_states[key]:
                        if key not in next_states and len(next_states)>=max_states:
                            raise ValueError('Joint envelope state budget exceeded')
                        next_states[key]=value
                        if return_witnesses:next_paths[key]=paths[(charge,squares,tail)]+(q,)
            if len(next_states)>max_states:
                raise ValueError('Joint envelope state budget exceeded')
            visited+=len(next_states);peak=max(peak,len(next_states));states=next_states
            paths=next_paths
        best=max((key for key in states if key[0]==0 and key[1]==2*D),key=states.__getitem__)
        maximum=states[best]
        lower=U*D-t*maximum/profile[D]
        rows.append({'doublons':D,'maximum_joint_penalty':str(maximum),'row_lower':str(lower),'target_margin':str(lower-gamma)})
        if return_witnesses:rows[-1]['worst_charge_pattern']=list(paths[best])
    return {'rows':rows,'minimum_lower':str(min(F(row['row_lower']) for row in rows)),
            'range':radius,'memory':memory,'local_table_entries':local_entries,'peak_states':peak,
            'visited_states':visited,'largest_local_window':largest_window,
            'local_window_spans_system':largest_window==sites,
            'count_profile':[str(value) for value in profile],
            'complete_neutral_charge_pattern_list_materialized':False,
            'scope':'Exact worst weighted Q row for a positive doublon-count profile times positive finite-range charge factors on the neutral open Hubbard chain. The profile vanishes at D=0. Dynamic programming merges prefixes and jointly sums local hopping penalties. Complexity remains exponential in interaction range; no general marginal-cone claim.'}


def replay(certificate, **budgets):
    receipt=joint_envelope(certificate, **budgets)
    if F(receipt['minimum_lower']) < F(certificate['target_lower']):
        raise ValueError('Joint charge envelope misses the requested complement bound')
    return dict(receipt,complement_lower=certificate['target_lower'],exact_accepted=True)


if __name__=='__main__':
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--verify',required=True);a=p.parse_args()
    print(json.dumps(replay(json.loads(Path(a.verify).read_text())),indent=2))
