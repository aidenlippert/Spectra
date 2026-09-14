"""Exact symmetric variational witness, with no configuration enumeration."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
import argparse
import json
from pathlib import Path

from experiments.marginal_symbolic import add, decode, mono, product, verify


def hopping_polynomial(modes, t, asymmetric=False):
    if type(modes) is not int or modes % 2 or modes < 4:
        raise ValueError("Model requires even M>=4")
    half = modes//2; t = F(t)
    n = [mono(((1, i), (0, i))) for i in range(modes)]
    h = add(*(product(n[i], n[j]) for group in (range(half), range(half, modes)) for i, j in combinations(group, 2)))
    edges = [(i, i+half, F(1)) for i in range(half)]
    if asymmetric:
        if modes != 6:
            raise ValueError("Asymmetric control is defined for M=6")
        edges = [(0, 3, F(1)), (1, 4, F(7, 10)), (2, 5, F(13, 10)), (0, 4, F(2, 5))]
    for i, j, weight in edges:
        h = add(h, mono(((1, i), (0, j)), -t*weight), mono(((1, j), (0, i)), -t*weight))
    return h


def upper_bound(modes, t, amplitudes):
    if type(modes) is not int or modes < 4 or modes % 2:
        raise ValueError("Even M>=4 required")
    half = modes//2; t = F(t)
    if len(amplitudes) != half+1 or any(type(a) is not int for a in amplitudes) or not any(amplitudes):
        raise ValueError("Expected nonzero integer symmetric amplitudes")
    norm = sum(comb(half,k)*a*a for k,a in enumerate(amplitudes))
    energy = sum(comb(half,k)*(k*(k-1)//2+(half-k)*(half-k-1)//2)*a*a
                 for k,a in enumerate(amplitudes))
    energy -= 2*t*sum(comb(half,k)*(half-k)*amplitudes[k]*amplitudes[k+1] for k in range(half))
    return {"norm": str(norm), "upper": str(F(energy,norm)), "upper_float": float(F(energy,norm))}


def verify_interval(certificate):
    modes = certificate["modes"]
    recipe = certificate["variational_upper"]
    if certificate["particles"] != modes//2:
        raise ValueError("Symmetric ansatz requires N=M/2")
    if not isinstance(recipe["t"], str):
        raise ValueError("t must be an exact rational string")
    t = F(recipe["t"])
    actual = decode(certificate["hamiltonian"], modes, 4)
    if actual != hopping_polynomial(modes,t):
        raise ValueError("Upper witness does not match the certified Hamiltonian")
    result = verify(certificate)
    result.update(upper_bound(modes,t,recipe["amplitudes"]))
    width = F(result["upper"])-F(result["lower"])
    if width < 0:
        raise ValueError("Inconsistent interval")
    result.update({"width": str(width), "width_float": float(width)})
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate",type=Path)
    args=parser.parse_args()
    print(json.dumps(verify_interval(json.loads(args.certificate.read_text())),indent=2))
