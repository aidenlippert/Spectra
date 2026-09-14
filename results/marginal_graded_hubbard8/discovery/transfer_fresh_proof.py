"""Bounded unchanged-proof transfer probe for the joint Reynolds certificate."""
import json
from fractions import Fraction as F
from pathlib import Path
from experiments.marginal_polynomial_metric import replay

ROOT = Path(__file__).parents[1]
SOURCE = ROOT / 'joint_reynolds_crossover/proof/certificate.json'
OUT = ROOT / 'fresh_transfer'


def perturb(data, epsilon, signs=(1, 1)):
    out = json.loads(json.dumps(data))
    pairs = []
    words = [tuple(tuple(x) for x in item['word']) for item in out['hamiltonian']]
    for word in words:
        if len(word) != 2 or word[0][0] != 1 or word[1][0] != 0:
            continue
        adj = ((1, word[1][1]), (0, word[0][1]))
        if adj in words and (word, adj) not in pairs and (adj, word) not in pairs:
            pairs.append((word, adj))
    if len(pairs) < 2:
        raise ValueError('Expected onsite and hopping probe terms')
    for sign, (word, adj) in zip(signs, pairs[:2]):
        for item in out['hamiltonian']:
            if tuple(tuple(x) for x in item['word']) in (word, adj):
                item['coefficient'] = str(F(item['coefficient']) + sign * epsilon)
    return out


def run():
    base = json.loads(SOURCE.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    accepted, rejected = [], []
    for den in (10**4, 10**5, 10**6):
        eps = F(1, den)
        for signs in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            name = f'eps_1e-{len(str(den))-1}_{"".join("p" if s>0 else "m" for s in signs)}'
            try:
                receipt = replay(perturb(base, eps, signs))
                record = {'epsilon': str(eps), 'signs': signs, 'receipt': receipt, 'scope': 'Exact replay with unchanged proof; Q complement bound, not ground energy.'}
                accepted.append(record); (OUT / (name + '.json')).write_text(json.dumps(record, indent=2) + '\n')
            except Exception as exc:
                record = {'epsilon': str(eps), 'signs': signs, 'rejection': str(exc), 'scope': 'Rejected unchanged-proof replay; Q complement bound, not ground energy.'}
                rejected.append(record); (OUT / ('rejected_' + name + '.json')).write_text(json.dumps(record, indent=2) + '\n')
    summary = {'source': str(SOURCE), 'accepted': accepted, 'rejected': rejected, 'scope': 'Bounded exact transfer probe; no solver and no claim of general transfer.'}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
