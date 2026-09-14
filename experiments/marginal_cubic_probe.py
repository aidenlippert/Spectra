"""Test all pure triples in addition to the full mixed-cubic dictionary."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from experiments.marginal_coefficient import dagger, export, solve_coefficients, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval
from experiments.marginal_compression import charge, compressed_dictionaries, hopping_components


def complete_cubic(h, modes):
    blocks = compressed_dictionaries(h,modes,'full',True)
    blocks = [b for b in blocks if not b['name'].startswith('local-triples')]
    groups = hopping_components(h,modes)
    membership = {i:g for g,group in enumerate(groups) for i in group}
    triples = [tuple((0,i) for i in indices) for indices in combinations(range(modes),3)]
    for name, words in [('all-triples-',triples),('all-triples+',[dagger(w) for w in triples])]:
        sectors = defaultdict(list)
        for w in words:
            sectors[charge(w,membership,len(groups))].append(w)
        blocks.extend({'name':f'{name}:{q}','words':ws} for q,ws in sorted(sectors.items()))
    return blocks


def run(modes=10,t=F(1,5)):
    h = hopping_polynomial(modes,t)
    blocks = complete_cubic(h,modes)
    solution = solve_coefficients(h,modes,modes//2,blocks)
    cert, receipt = export(h,modes,modes//2,blocks,solution)
    cert['variational_upper'] = symmetric_upper(modes,t)
    receipt.update(verify_interval(cert))
    out = Path(__file__).resolve().parents[1]/'results/marginal_sector_reference'
    out.mkdir(exist_ok=True)
    (out/'complete_cubic.json').write_text(json.dumps(cert)+'\n')
    (out/'complete_cubic_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__ == '__main__': run()
