"""Discover fermionic SOS certificates directly in CAR coefficient space.

No Fock basis, sector matrices, old certificates, or many-body eigensolver are
used by discovery. Numerical PSD optimization proposes factors; the existing
exact symbolic checker alone determines the accepted lower bound.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import argparse
import json
import time

import numpy as np
import scipy.sparse as sp
import cvxpy as cp

from experiments.marginal_symbolic import (
    add, canonical, encode, hermitian, mono, multiplier_basis, number_shift,
    product, scale, verify, word_product,
)
from experiments.marginal_collective import hopping_polynomial, upper_bound


def dagger(word):
    return tuple((1-c, i) for c, i in reversed(word))


def dictionaries(modes, family):
    if modes % 2 or modes < 4:
        raise ValueError("Test dictionaries require even M>=4")
    if family not in ("quadratic", "local", "mixed"):
        raise ValueError("Unknown operator family")
    linear = [((0, i),) for i in range(modes)]
    pairs = [((0, j), (0, i)) for i, j in combinations(range(modes), 2)]
    number = [((1, i), (0, j)) for i in range(modes) for j in range(modes)]
    blocks = [{"name": name, "words": words} for name, words in
              (("linear-", linear), ("linear+", [dagger(w) for w in linear]),
               ("pair-", pairs), ("pair+", [dagger(w) for w in pairs]),
               ("particle-hole", number))]
    if family != "quadratic":
        half = modes // 2
        groups = (range(half), range(half, modes))
        triple = [tuple((0, i) for i in indices) for group in groups for indices in combinations(group, 3)]
        if triple:
            blocks += [{"name": "local-triples-", "words": triple},
                       {"name": "local-triples+", "words": [dagger(w) for w in triple]}]
        mixed = list(linear)
        if family == "mixed":
            mixed += [((1, k), (0, j), (0, i)) for i, j in combinations(range(modes), 2) for k in range(modes)]
        else:
            mixed += [((1, k), (0, j), (0, i)) for group in groups for i, j in combinations(group, 2) for k in group]
        blocks += [{"name": f"{family}-", "words": mixed},
                   {"name": f"{family}+", "words": [dagger(w) for w in mixed]}]
    return blocks


def coefficient_rows(modes):
    # A real Hermitian balanced polynomial has symmetric coefficient pairs.
    # Match one from each pair; exact replay later checks both partners.
    rows = []
    for k in range(4):
        sets = list(combinations(range(modes), k))
        for i, left in enumerate(sets):
            for right in sets[i:]:
                rows.append(tuple((1, p) for p in left) + tuple((0, p) for p in right))
    return rows


def gram_map(words, lookup):
    n = len(words)
    row_indices, cols, values = [], [], []
    for i, left in enumerate(words):
        for j, right in enumerate(words):
            for w, coefficient in word_product(dagger(left), right):
                if w in lookup:
                    row_indices.append(lookup[w])
                    cols.append(i*n+j)
                    values.append(float(coefficient))
    return sp.csr_matrix((values, (row_indices, cols)), shape=(len(lookup), n*n))


def solve_coefficients(h, modes, particles, blocks, use_number_ideal=True):
    if type(modes) is not int or modes < 1 or type(particles) is not int or not 0 <= particles <= modes:
        raise ValueError("Invalid sector")
    h = canonical(h)
    if not hermitian(h) or any(len(w)>4 or sum(2*c-1 for c, _ in w) for w in h):
        raise ValueError("Expected a real Hermitian two-body Hamiltonian")
    if any(i not in range(modes) for w in h for _, i in w):
        raise ValueError("Hamiltonian index outside declared mode count")
    start = time.monotonic()
    rows = coefficient_rows(modes)
    lookup = {w: i for i, w in enumerate(rows)}
    basis = multiplier_basis(modes) if use_number_ideal else []
    shift = number_shift(modes, particles)
    columns = [product(shift, q) for q in basis]
    a = sp.csr_matrix(np.array([[float(col.get(w, 0)) for col in columns] for w in rows])) if basis else None
    b = cp.Variable()
    unit = np.zeros(len(rows)); unit[lookup[()]] = 1
    expression = b*unit
    xvar = cp.Variable(len(basis)) if basis else None
    if basis:
        expression = expression + a @ xvar
    variables, shapes, nonzeros = [], [], 0
    for block in blocks:
        words = block["words"]
        if not words or len({sum(2*c-1 for c, _ in w) for w in words}) != 1:
            raise ValueError("Gram dictionary must have a single charge")
        if any(len(w)>3 or any(c not in (0, 1) or i not in range(modes) for c, i in w) for w in words):
            raise ValueError("Invalid dictionary word")
        g = gram_map(words, lookup)
        q = cp.Variable((len(words), len(words)), PSD=True)
        expression = expression + g @ cp.reshape(q, (len(words)**2,), order="C")
        variables.append(q); shapes.append(len(words)); nonzeros += g.nnz
    constraint = expression == np.array([float(h.get(w, 0)) for w in rows])
    problem = cp.Problem(cp.Maximize(b), [constraint])
    build_seconds = time.monotonic()-start
    start = time.monotonic()
    problem.solve(solver="SCS", eps=1e-8, max_iters=30000)
    solve_seconds = time.monotonic()-start
    if b.value is None or any(q.value is None for q in variables) or (basis and xvar.value is None):
        raise RuntimeError(f"No coefficient certificate proposal: {problem.status}")
    return {"b": float(b.value), "grams": [q.value for q in variables],
            "x": np.asarray(xvar.value) if basis else np.array([]), "basis": basis,
            "status": problem.status, "coefficient_rows": len(rows),
            "gram_dimensions": shapes, "coefficient_map_nonzeros": nonzeros,
            "multiplier_basis_dimension": len(basis), "build_seconds": build_seconds,
            "solve_seconds": solve_seconds}


def export(h, modes, particles, blocks, solution, denominator=10**7, operator_degree=3):
    factors = []
    for idx, (block, q) in enumerate(zip(blocks, solution["grams"])):
        direct = solution.get("factor_matrices", [None]*len(blocks))[idx]
        if direct is None:
            values, vectors = np.linalg.eigh((q+q.T)/2)
            l = np.sqrt(np.maximum(values, 0))[:, None]*vectors.T
        else:
            l = np.asarray(direct)
        integers = [[int(round(value*denominator)) for value in row] for row in l]
        integers = [row for row in integers if any(row)]
        factors.append({"name": block["name"], "words": block["words"], "factor": integers})
    x = add(*(scale(q, F(int(round(value*10**10)), 10**10)) for q, value in zip(solution["basis"], solution["x"])))
    certificate = {"modes": modes, "particles": particles, "hamiltonian": encode(h),
                   "number_multiplier": encode(x), "b": str(F(int(round(solution["b"]*10**12)), 10**12)),
                   "denominator": denominator, "blocks": factors}
    if operator_degree != 3:
        certificate["operator_degree"] = operator_degree
    receipt = verify(certificate)
    receipt.update({k: solution[k] for k in ("status", "coefficient_rows", "gram_dimensions",
                    "coefficient_map_nonzeros", "multiplier_basis_dimension", "build_seconds", "solve_seconds")})
    receipt["proposed_b"] = solution["b"]
    return certificate, receipt


def hopping_model(modes, t, asymmetric=False):
    return hopping_polynomial(modes,t,asymmetric)


def symmetric_upper(modes, t, denominator=10**8):
    """Small symmetric variational ansatz, not full-sector diagonalization.

One particle per matched flavor pair. a[k] is a common amplitude for every
flavor-ordered configuration with k particles on the left. Combinatorial
counts give its norm and energy without enumerating those configurations.
"""
    half = modes//2; t = F(t)
    diagonal = [k*(k-1)//2 + (half-k)*(half-k-1)//2 for k in range(half+1)]
    mat = np.diag(np.array(diagonal, dtype=float))
    for k in range(half):
        mat[k,k+1] = mat[k+1,k] = -float(t)*np.sqrt((k+1)*(half-k))
    _, vectors = np.linalg.eigh(mat)
    amplitudes = [int(round(vectors[k,0]/np.sqrt(comb(half,k))*denominator)) for k in range(half+1)]
    return {"ansatz": "one particle per matched flavor; permutation-symmetric amplitudes",
            "amplitudes": amplitudes, "t": str(t), "variational_dimension": half+1,
            **upper_bound(modes,t,amplitudes)}


def run_case(modes, t, family, asymmetric=False, use_number_ideal=True):
    h = hopping_model(modes,t,asymmetric)
    blocks = dictionaries(modes,family)
    solution = solve_coefficients(h,modes,modes//2,blocks,use_number_ideal)
    certificate, receipt = export(h,modes,modes//2,blocks,solution)
    receipt.update({"modes": modes, "t": str(t), "family": family, "asymmetric": asymmetric,
                    "use_number_ideal": use_number_ideal})
    if not asymmetric:
        upper = symmetric_upper(modes,t)
        certificate["variational_upper"] = upper
        receipt.update({"upper": upper["upper"], "width": str(F(upper["upper"])-F(receipt["lower"])),
                        "width_float": float(F(upper["upper"])-F(receipt["lower"]))})
    return certificate,receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes",type=int,default=6)
    parser.add_argument("--t",default="1/5")
    parser.add_argument("--family",choices=("quadratic","local","mixed"),default="mixed")
    parser.add_argument("--asymmetric",action="store_true")
    parser.add_argument("--no-number-ideal",action="store_true")
    args=parser.parse_args()
    certificate,receipt=run_case(args.modes,F(args.t),args.family,args.asymmetric,not args.no_number_ideal)
    out=Path(__file__).resolve().parents[1]/"results/marginal_coefficient";out.mkdir(exist_ok=True)
    name=f"m{args.modes}_t{str(F(args.t)).replace('/','_')}_{args.family}"+('_asymmetric' if args.asymmetric else '')+('_no_ideal' if args.no_number_ideal else '')
    (out/f"{name}.json").write_text(json.dumps(certificate)+'\n')
    receipt["certificate_file"]=f"{name}.json"
    (out/f"{name}_receipt.json").write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=="__main__":
    main()
