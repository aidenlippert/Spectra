"""Exact three-spectator consistency tests beyond the two-spectator family."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_local_hubbard_block import _reflection, _sector
from experiments.marginal_charged_projectors import _particlehole, _spinflip, charged_vectors
from experiments.marginal_polynomial_sos import integer_psd
from experiments.marginal_spectator_hopping import LABELS as OLD_LABELS, actions as old_actions
from experiments.marginal_two_spectator_hopping import LABELS as TWO_LABELS, actions as two_actions


def reflected(label):
    i,j,k,l,m,p,q,r = label
    return (4-j,4-i,4-m,4-l,4-k,r,q,p)


LABELS = tuple(label for i,j in combinations(range(5),2)
    for k,l,m in [sorted(set(range(5))-{i,j})]
    for p,q,r in product((1,2),repeat=3)
    if (p+q+r+j-i)%2==1
    for label in [(i,j,k,l,m,p,q,r)] if label < reflected(label))


def components(label):
    mirror=reflected(label)
    shift=lambda t:tuple(v+1 for v in t[:5])+t[5:]
    return [(*label,1),(*mirror,-1),(*shift(label),-1),(*shift(mirror),1)]


def operator(label):
    terms=components(label); result=[]
    for state in range(4096):
        charges=[((state>>(2*i))&3).bit_count()-1 for i in range(6)]
        row={}
        for i,j,k,l,m,p,q,r,sign in terms:
            weight=sign*charges[k]**p*charges[l]**q*charges[m]**r
            for a,b in ((i,j),(j,i)):
                for spin in (0,1):
                    target=apply_word(((1,2*a+spin),(0,2*b+spin)),state)
                    if target is not None:
                        t,phase=target;row[t]=row.get(t,0)+weight*phase
        result.append({t:a for t,a in row.items() if a})
    for s,row in enumerate(result):
        for t,a in row.items():
            if _sector(s,6)!=_sector(t,6) or (s^t).bit_count()!=2 or result[t].get(s,0)!=a:
                raise ValueError('Physical hopping action failed')
        for fn in (lambda s:_reflection(s,6),_particlehole,_spinflip):
            us,phase=fn(s)
            right={fn(t)[0]:fn(t)[1]*a for t,a in row.items()}
            if right!={t:phase*a for t,a in result[us].items()}:
                raise ValueError('Three-spectator symmetry failed')
    for sites in (8,9,10,17):
        totals={}
        for offset in range(sites):
            for i,j,k,l,m,p,q,r,sign in terms:
                pair=tuple(sorted(((i+offset)%sites,(j+offset)%sites)))
                spectators=tuple(sorted((((k+offset)%sites,p),((l+offset)%sites,q),((m+offset)%sites,r))))
                key=pair+spectators;totals[key]=totals.get(key,0)+sign
        if any(totals.values()):raise ValueError('Periodic telescope failed')
    return result


def constant_operator(i,j):
    data=[]
    for s in range(4096):
        row={}
        for a,b,weight in [(i,j,1),(4-j,4-i,-1),(i+1,j+1,-1),(5-j,5-i,1)]:
            for left,right in ((a,b),(b,a)):
                for spin in (0,1):
                    target=apply_word(((1,2*left+spin),(0,2*right+spin)),s)
                    if target is not None:
                        t,phase=target;row[t]=row.get(t,0)+weight*phase
        data.append({t:a for t,a in row.items() if a})
    return data


def hopping_norm_bound():
    action=[]
    for state in range(16):
        row={}
        for i,j in ((0,1),(1,0)):
            for spin in (0,1):
                target=apply_word(((1,2*i+spin),(0,2*j+spin)),state)
                if target is not None:
                    s,phase=target;row[s]=row.get(s,0)+phase
        action.append({s:a for s,a in row.items() if a})
    if any(action[t].get(s,0)!=a for s,row in enumerate(action) for t,a in row.items()):
        raise ValueError('Two-site hopping is not Hermitian')
    bound=max(sum(abs(a) for a in row.values()) for row in action)
    if bound!=2:raise ValueError('Two-site hopping row-sum bound changed')
    return bound


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    folder=parser.parse_args().directory.resolve()
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json'
    c=json.loads(cp.read_text());r=json.loads(rp.read_text())
    if not r.get('accepted') or not r.get('two_spectator_hopping'):
        raise ValueError('Accepted two-spectator family required')
    for name,expected in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Stale accepted source: '+name)
    if len(LABELS)!=18:raise ValueError('Expected full18-direction three-spectator basis')
    base_norm=hopping_norm_bound()
    data=[operator(label) for label in LABELS]
    old=[constant_operator(*pair) for pair in ((0,1),(0,3),(1,2))]+[
        old_actions({label:1}) for label in OLD_LABELS]+[two_actions({label:1}) for label in TWO_LABELS]
    # Restrict to off-diagonal two-bit changes with total particle number<=4.
    # Every fixed projector vanishes there; diagonal, spin and pair operators
    # have no such entries. Retain all three unconditioned hopping directions,
    # including the two variable nearest-neighbor profiles and thirty two-spectator terms.
    charged,_=charged_vectors(c['charged_vector'])
    sectors={int(s).bit_count() for s,a in c['half_vector'].items() if a}|{
        s.bit_count() for v in charged for s,a in v.items() if a}
    if sectors!={5,6,7}:raise ValueError('Fixed projector sectors changed')
    columns=[{(s,t):a for s,row in enumerate(action) if s.bit_count()<=4
              for t,a in row.items() if t>s} for action in old+data]
    gram=[[sum(a*right.get(key,0) for key,a in left.items()) for right in columns] for left in columns]
    old_rank=integer_psd([row[:len(old)] for row in gram[:len(old)]])
    total_rank=integer_psd(gram)
    moments=[]
    for label,action in zip(LABELS,data):
        value=F(0)
        for item in c['mixture']:
            vector={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in vector.values())
            value+=F(item['weight'])*F(sum(a*b*vector.get(t,0) for s,a in vector.items() for t,b in action[s].items()),norm)
        moments.append({'label':list(label),'exact_moment':str(value),'moment_float':float(value),
                        'norm_upper_bound':8,'normalized_violation_lower':float(abs(value)/8),'violated':bool(value)})
    files={Path(__file__).resolve(),cp,rp}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'directions':18,'violated_directions':sum(m['violated'] for m in moments),
            'two_site_hopping_row_sum_bound':base_norm,
            'moments':moments,'old_rank':old_rank,'combined_rank':total_rank,
            'new_independent_directions':total_rank['rank']-old_rank['rank'],'restricted_integer_gram':gram,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'Exact necessary stationary consistency tests for18 PH-even reflected three-spectator charge-hopping telescopes. Full-Fock CAR, Hermiticity, reflection/PH/spin-flip symmetry and periodic cancellation checked. Each term is a product of three disjoint charge powers (norm<=1) and spin-summed hopping (norm2); triangle norm bound8. Restricted integer Gram measures independence modulo all47 preceding hopping directions; other preceding operators have no entries on this restriction. Nonzero symmetry-invariant moments exclude a stationary extension of this particular averaged local mixture. No stronger energy bound, enlarged-family optimum or general representability claim.'}
    (folder/'three_spectator_overlap.json').write_text(json.dumps(result,indent=2)+'\n')
    largest=sorted(moments,key=lambda m:abs(m['moment_float']),reverse=True)[:4]
    print(json.dumps({'accepted':True,'directions':18,'violated':result['violated_directions'],
                      'old_rank':old_rank['rank'],'combined_rank':total_rank['rank'],
                      'largest':[{'label':m['label'],'moment':m['moment_float']} for m in largest]}))


if __name__=='__main__':main()
