"""Exact classical Markov extension of the binary projection of an accepted dual."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parents[3]


def extend(probabilities):
    if len(probabilities)!=64 or any(p<0 for p in probabilities) or sum(probabilities)!=1:
        raise ValueError('Normalized nonnegative six-bit distribution required')
    prefix=[sum(probabilities[s] for s in range(64) if s&31==x) for x in range(32)]
    suffix=[sum(probabilities[s] for s in range(64) if s>>1==x) for x in range(32)]
    if prefix!=suffix:raise ValueError('Five-bit marginals disagree')
    transition=[[F(0)]*32 for _ in range(32)]
    for word,p in enumerate(probabilities):
        if prefix[word&31]:transition[word&31][word>>1]+=p/prefix[word&31]
    for x,p in enumerate(prefix):
        if not p:transition[x][x>>1]=F(1)
    if any(sum(row)!=1 or any(v<0 for v in row) for row in transition):raise ValueError('Invalid stochastic matrix')
    if any(transition[x][y] and y&15!=x>>1 for x in range(32) for y in range(32)):raise ValueError('Markov transition violates the shift')
    stationary=[sum(prefix[x]*transition[x][y] for x in range(32)) for y in range(32)]
    if stationary!=prefix:raise ValueError('Stationary measure fails')
    if any(prefix[word&31]*transition[word&31][word>>1]!=p for word,p in enumerate(probabilities)):raise ValueError('Six-bit distribution is not reconstructed')
    return prefix,transition


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve()
    # Independent known models, plus refusal of inconsistent prefix/suffix data.
    stationary,transition=extend([F(1,64)]*64)
    assert stationary==[F(1,32)]*32 and all(sorted(v for v in row if v)==[F(1,2)]*2 for row in transition)
    alternating=[F(0)]*64;alternating[21]=alternating[42]=F(1,2)
    stationary,transition=extend(alternating)
    assert stationary[10]==stationary[21]==F(1,2) and transition[10][21]==transition[21][10]==1
    inconsistent=[F(0)]*64;inconsistent[1]=F(1)
    try:extend(inconsistent)
    except ValueError:pass
    else:raise AssertionError('Inconsistent distribution accepted')
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json'
    receipt=json.loads(rp.read_text())
    if not receipt['accepted'] or not receipt.get('full_charge_indicators') or receipt['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted full indicator family required')
    certificate=json.loads(cp.read_text());probabilities=[F(0)]*64
    for item in certificate['mixture']:
        vector={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in vector.values());counts=[0]*64
        for s,a in vector.items():
            mask=sum(int(((s>>(2*i))&3) in (0,3))<<i for i in range(6));counts[mask]+=a*a
        for mask,total in enumerate(counts):probabilities[mask]+=F(item['weight'])*F(total,norm)
    stationary,transition=extend(probabilities)
    result={'accepted':True,'classical_only':True,'known_model_checks_passed':3,'six_bit_probabilities':list(map(str,probabilities)),
        'stationary_five_bit_probabilities':list(map(str,stationary)),
        'transitions':[[x,y,str(v)] for x,row in enumerate(transition) for y,v in enumerate(row) if v],
        'positive_stationary_states':sum(bool(p) for p in stationary),'positive_six_bit_patterns':sum(bool(p) for p in probabilities),
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),cp,rp]},
        'scope':'Exact stationary order-five classical Markov extension of the six-site binary empty/double indicator distribution of this local quantum mixture. All 32 prefix/suffix equalities, stochastic rows, shift compatibility, stationarity and all 64 reconstructed probabilities are checked. Infinite stationary classical sequences follow by iterating this finite stochastic kernel. Signed charges, spins, fermionic coherences and extension of the original quantum density matrix are not certified.'}
    (folder/'indicator_markov_extension.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','six_bit_probabilities','stationary_five_bit_probabilities','transitions')}))


if __name__=='__main__':main()
