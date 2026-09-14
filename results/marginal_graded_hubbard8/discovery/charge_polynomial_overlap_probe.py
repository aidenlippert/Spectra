"""Exact PH-even charge-polynomial consistency separators in a certified dual."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]


def monomial(q, exponents):
    return prod(x**power for x, power in zip(q, exponents))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    folder = parser.parse_args().directory.resolve()
    cp = folder/'range_two_family_limit_certificate.json'
    rp = folder/'range_two_family_limit_replay.json'
    receipt = json.loads(rp.read_text())
    if not receipt['accepted'] or not receipt.get('full_quadratic_charge') or receipt['source_sha256'][str(cp.relative_to(ROOT))] != hashlib.sha256(cp.read_bytes()).hexdigest():
        raise ValueError('Accepted full quadratic family mixture required')
    c = json.loads(cp.read_text())
    states = []
    for item in c['mixture']:
        vector = {int(s): a for s, a in item['vector'].items()}
        states.append((F(item['weight']), vector, sum(a*a for a in vector.values())))
    qs = [[((s >> (2*i)) & 3).bit_count()-1 for i in range(6)] for s in range(4096)]
    directions = []
    for exponents in product(range(3), repeat=5):
        reflected = exponents[::-1]
        if sum(exponents) % 2 or exponents >= reflected:
            continue
        # q in {-1,0,1}, so powers 0,1,2 span every function of one charge.
        ys = [monomial(q[:5], exponents)-monomial(q[:5], reflected) for q in qs[:1024]]
        ts = [ys[s & 1023]-ys[s >> 2] for s in range(4096)]
        # Independent direct six-site expansion and PH/reflection identities.
        if any(ts[s] != monomial(q[:5],exponents)-monomial(q[:5],reflected)-monomial(q[1:],exponents)+monomial(q[1:],reflected) for s,q in enumerate(qs)):
            raise ValueError('Physical shift disagrees with polynomial expansion')
        for s in range(1024):
            rev = sum(((s >> (2*i)) & 3) << (2*(4-i)) for i in range(5))
            if ys[rev] != -ys[s] or ys[s ^ 1023] != ys[s]:
                raise ValueError('Reflection or particle-hole symmetry fails')
        moment = sum((w*F(sum(a*a*ts[s] for s,a in vector.items()),norm) for w,vector,norm in states),F(0))
        if sum(exponents) == 2 and moment:
            raise ValueError('Previously certified quadratic moment is nonzero')
        norm = max(map(abs,ts))
        directions.append({'exponents':exponents,'reflected_exponents':reflected,'degree':sum(exponents),
            'exact_moment':str(moment),'moment_float':float(moment),'local_operator_norm':norm,
            'normalized_violation':float(abs(moment)/norm) if norm else 0,
            'five_site_nonzero_entries':sum(bool(v) for v in ys)})
    if len(directions) != 52:
        raise ValueError('Expected the complete 52-element PH-even reflection-odd charge basis')
    directions.sort(key=lambda row:(-row['normalized_violation'],row['exponents']))
    result={'accepted':True,'directions':directions,'basis_size':len(directions),
        'violated_directions':sum(F(row['exact_moment'])!=0 for row in directions),
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),cp,rp]},
        'scope':'Exact necessary translation-consistency moments on an accepted local PSD mixture. Complete PH-even reflection-odd basis of functions of five site charges, with each single-site power at most two. No optimization or energy improvement is certified by this probe. Does not include spin-resolved or off-diagonal operators.'}
    (folder/'charge_polynomial_overlap.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'basis_size':len(directions),'violated_directions':result['violated_directions'],'top':[ {k:v for k,v in row.items() if k!='exact_moment'} for row in directions[:6]]}),flush=True)


if __name__ == '__main__':
    main()
