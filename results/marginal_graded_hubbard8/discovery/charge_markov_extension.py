"""Exact stationary classical extension of the symmetry-averaged signed charges."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib
ROOT=Path(__file__).resolve().parents[3]


def reverse(word,digits):
    out=0
    for _ in range(digits):out=3*out+word%3;word//=3
    return out


def extend(probabilities):
    if len(probabilities)!=729 or sum(probabilities)!=1 or any(p<0 for p in probabilities):raise ValueError('Normalized signed-charge distribution required')
    pi=[F(0)]*243;suffix=[F(0)]*243
    for word,p in enumerate(probabilities):pi[word%243]+=p;suffix[word//3]+=p
    if pi!=suffix:raise ValueError('Five-charge prefix/suffix marginals differ')
    rows=[{} for _ in range(243)]
    for word,p in enumerate(probabilities):
        if p:rows[word%243][word//3]=p/pi[word%243]
    for x,p in enumerate(pi):
        if not p:rows[x][x//3]=F(1)
    stationary=[F(0)]*243
    for x,row in enumerate(rows):
        if sum(row.values())!=1 or any(v<0 for v in row.values()):raise ValueError('Invalid stochastic row')
        for y,v in row.items():
            if y%81!=x//3:raise ValueError('Transition is not a shift')
            stationary[y]+=pi[x]*v
    if stationary!=pi:raise ValueError('Stationarity fails')
    if any(pi[word%243]*rows[word%243].get(word//3,F(0))!=p for word,p in enumerate(probabilities)):raise ValueError('Six-charge law is not recovered')
    return pi,rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);folder=parser.parse_args().directory.resolve()
    pi,rows=extend([F(1,729)]*729);assert pi==[F(1,243)]*243 and all(list(row.values())==[F(1,3)]*3 for row in rows)
    bad=[F(0)]*729;bad[1]=F(1)
    try:extend(bad)
    except ValueError:pass
    else:raise AssertionError('Inconsistent law accepted')
    cp=folder/'range_two_family_limit_certificate.json';rp=folder/'range_two_family_limit_replay.json';receipt=json.loads(rp.read_text())
    if not receipt['accepted'] or not receipt.get('full_signed_charge') or receipt['source_sha256'][str(cp.relative_to(ROOT))]!=hashlib.sha256(cp.read_bytes()).hexdigest():raise ValueError('Accepted full signed-charge family required')
    c=json.loads(cp.read_text());raw=[F(0)]*729
    for item in c['mixture']:
        vector={int(s):a for s,a in item['vector'].items()};norm=sum(a*a for a in vector.values());counts=[0]*729
        for s,a in vector.items():counts[sum(((s>>(2*i))&3).bit_count()*3**i for i in range(6))]+=a*a
        for word,n in enumerate(counts):raw[word]+=F(item['weight'])*F(n,norm)
    symmetric=[(raw[word]+raw[728-word]+raw[reverse(word,6)]+raw[728-reverse(word,6)])/4 for word in range(729)]
    pi,rows=extend(symmetric)
    raw_prefix=[F(0)]*243;raw_suffix=[F(0)]*243
    for word,p in enumerate(raw):raw_prefix[word%243]+=p;raw_suffix[word//3]+=p
    result={'accepted':True,'classical_only':True,'symmetrization':'Equal mixture of identity, charge sign reversal, spatial reflection and their composition on the six-charge law.','known_model_checks_passed':2,'raw_prefix_suffix_l1_difference':str(sum(abs(a-b) for a,b in zip(raw_prefix,raw_suffix))),'six_charge_probabilities':list(map(str,symmetric)),'stationary_five_charge_probabilities':list(map(str,pi)),'transitions':[[x,y,str(v)] for x,row in enumerate(rows) for y,v in row.items()],'positive_stationary_states':sum(bool(p) for p in pi),'positive_six_charge_patterns':sum(bool(p) for p in symmetric),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),cp,rp]},'scope':'Exact order-five stationary classical Markov extension of the PH/reflection-averaged signed-charge projection. All243prefix/suffix equalities, stochastic rows, shift compatibility, stationarity and729reconstructed probabilities checked. This does not certify fixed global particle number, spins, coherences or extension of the quantum density matrix.'}
    (folder/'charge_markov_extension.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','six_charge_probabilities','stationary_five_charge_probabilities','transitions')}))


if __name__=='__main__':main()
