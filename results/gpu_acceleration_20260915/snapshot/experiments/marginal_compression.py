"""Compress cubic dictionaries and split Grams by conserved hopping charges."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import json
from pathlib import Path

from experiments.marginal_coefficient import dictionaries, export, solve_coefficients, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval


def hopping_components(h, modes):
    """Components of one-body hopping; reject other terms that violate them."""
    neighbors = [set() for _ in range(modes)]
    for word, coefficient in h.items():
        if coefficient and len(word) == 2 and word[0][0] == 1 and word[1][0] == 0:
            i, j = word[0][1], word[1][1]
            neighbors[i].add(j); neighbors[j].add(i)
    groups = []
    remaining = set(range(modes))
    while remaining:
        pending = [min(remaining)]; group = set()
        while pending:
            i = pending.pop()
            if i in group:
                continue
            group.add(i); pending.extend(neighbors[i]-group)
        remaining -= group; groups.append(sorted(group))
    membership = {i: g for g, group in enumerate(groups) for i in group}
    for word, coefficient in h.items():
        if coefficient and any(charge(word, membership, len(groups))):
            raise ValueError("Hamiltonian violates the proposed component charges")
    return groups


def charge(word, membership, count):
    result = [0]*count
    for creation, mode in word:
        result[membership[mode]] += 2*creation-1
    return tuple(result)


def compressed_dictionaries(h, modes, subset="full", split=False):
    if subset not in ("full", "density", "edge"):
        raise ValueError("Unknown cubic subset")
    blocks = dictionaries(modes, "mixed")
    edges = {frozenset((w[0][1], w[1][1])) for w, c in h.items()
             if c and len(w) == 2 and w[0][1] != w[1][1]}
    if subset != "full":
        # Select in the charge-minus block and derive its adjoint partner by
        # matching support. k repeated among i,j gives a density-dressed linear.
        for block in blocks[-2:]:
            selected = []
            for w in block["words"]:
                if len(w) == 1:
                    selected.append(w); continue
                k = w[0][1] if sum(2*c-1 for c, _ in w) == -1 else w[-1][1]
                pair = [i for c, i in w if c == (0 if sum(2*c-1 for c, _ in w) == -1 else 1)]
                if k in pair or (subset == "edge" and any(frozenset((k,i)) in edges for i in pair)):
                    selected.append(w)
            block["words"] = selected
    if split:
        groups = hopping_components(h, modes)
        membership = {i:g for g, group in enumerate(groups) for i in group}
        result = []
        for block in blocks:
            sectors = defaultdict(list)
            for w in block["words"]:
                sectors[charge(w, membership, len(groups))].append(w)
            result.extend({"name": f"{block['name']}:{q}", "words": words} for q, words in sorted(sectors.items()))
        blocks = result
    return blocks


def run(modes, t, subset, split, asymmetric=False):
    h = hopping_polynomial(modes, t, asymmetric)
    blocks = compressed_dictionaries(h, modes, subset, split)
    solution = solve_coefficients(h, modes, modes//2, blocks)
    certificate, receipt = export(h, modes, modes//2, blocks, solution)
    if not asymmetric:
        certificate["variational_upper"] = symmetric_upper(modes,t)
        receipt.update(verify_interval(certificate))
    receipt.update({"modes": modes, "t": str(t), "subset": subset, "split": split,
                    "asymmetric": asymmetric, "gram_scalar_variables": sum(len(b["words"])*(len(b["words"])+1)//2 for b in blocks)})
    if split:
        receipt["conserved_components"] = hopping_components(h,modes)
    return certificate, receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, default=8)
    parser.add_argument("--t", default="1/5")
    parser.add_argument("--subset", choices=("full", "density", "edge"), default="full")
    parser.add_argument("--split", action="store_true")
    parser.add_argument("--asymmetric", action="store_true")
    args = parser.parse_args()
    cert, receipt = run(args.modes,F(args.t),args.subset,args.split,args.asymmetric)
    out = Path(__file__).resolve().parents[1]/"results/marginal_compression"
    out.mkdir(exist_ok=True)
    name = f"m{args.modes}_t{str(F(args.t)).replace('/','_')}_{args.subset}"+('_split' if args.split else '')+('_asymmetric' if args.asymmetric else '')
    (out/f"{name}.json").write_text(json.dumps(cert)+'\n')
    receipt["certificate_file"] = f"{name}.json"
    (out/f"{name}_receipt.json").write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__ == "__main__":
    main()
