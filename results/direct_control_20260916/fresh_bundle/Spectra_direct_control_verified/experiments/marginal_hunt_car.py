"""Exact small CAR polynomial diagnostic for cubic anticommutators."""
from fractions import Fraction
import json
from pathlib import Path

Word = tuple[tuple[int, int], ...]  # (creation: 0/1, mode)

def add(*ps):
    out = {}
    for p in ps:
        for w, c in p.items(): out[w] = out.get(w, Fraction(0)) + c
    return {w:c for w,c in out.items() if c}

def scale(p, x): return {w:c*x for w,c in p.items() if c*x}

def mul(p, q):
    out = {}
    for w, a in p.items():
        for v, b in q.items():
            z = {(w+v): a*b}
            # recurse on first inversion (annihilator left of creator)
            changed = True
            while changed:
                changed = False
                nz = {}
                for x, c in z.items():
                    for i in range(len(x)-1):
                        if x[i][0] == 0 and x[i+1][0] == 1:
                            changed = True
                            # a_i a_j^ = delta_ij - a_j^ a_i
                            if x[i][1] == x[i+1][1]:
                                y = x[:i] + x[i+2:]
                                nz[y] = nz.get(y, 0) + c
                            y = x[:i] + (x[i+1], x[i]) + x[i+2:]
                            nz[y] = nz.get(y, 0) - c
                            break
                    else:
                        nz[x] = nz.get(x, 0) + c
                z = {x:c for x,c in nz.items() if c}
            # canonical order and Pauli nilpotence
            for x,c in z.items():
                key = lambda q: (0 if q[0] else 1, q[1])
                inv = sum(key(x[i]) > key(x[j]) for i in range(len(x)) for j in range(i+1,len(x)))
                y = tuple(sorted(x, key=key))
                # Pauli nilpotence applies to equal operators after sorting,
                # even when equal modes were separated in the input word.
                if any(y[i] == y[i+1] for i in range(len(y)-1)): continue
                out[y] = out.get(y, 0) + c*((-1)**inv)
    return {w:c for w,c in out.items() if c}

def adj(p):
    return {tuple((1-int(cr), m) for cr,m in reversed(w)): c for w,c in p.items()}

def mono(items, c=1): return {tuple(items): Fraction(c)}

def run(kind):
    # Distinct labels prevent accidental nilpotence; rational sparse coefficients.
    if kind == "T1_charge_-3":
        b = add(mono(((0,0),(0,1),(0,2)), 2), mono(((0,0),(0,2),(0,3)), -1))
    elif kind == "T2_charge_-1": # one creator, two annihilators
        b = add(mono(((1,0),(0,1),(0,2)), 2), mono(((1,1),(0,2),(0,3)), -1))
    else:
        raise ValueError("Unknown certificate family")
    anti = add(mul(b, adj(b)), mul(adj(b), b))
    bb = mul(adj(b), b)
    degrees = sorted(set(len(w) for w in anti))
    high = {w:c for w,c in anti.items() if len(w) >= 6}
    return {"kind":kind, "anticommutator_degrees":degrees,
            "degree6_terms":len(high), "remainder_max_degree":max(degrees, default=0),
            "single_square_has_degree6":any(len(w)>=6 for w in bb),
            "coefficients":{str(w):str(c) for w,c in anti.items()}}

if __name__ == "__main__":
    result = {k:run(k) for k in ("T1_charge_-3","T2_charge_-1")}
    result_path = Path(__file__).resolve().parents[1] / "results/marginal_hunt_car.json"
    result_path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:{q:v for q,v in r.items() if q!="coefficients"} for k,r in result.items()}, indent=2))
