"""Bounded exact Pauli closure scaling experiment."""
from __future__ import annotations
import json, time
from collections import deque
from pathlib import Path
from fractions import Fraction
from .pauli import commutator_i, generator_matrix

def label(n, site, p):
    return "".join(p if i == site else "I" for i in range(n))

def quadratic_generators(n):
    return [{label(n, i, "Z"): Fraction(1)} for i in range(n)] + [{xx_label(n, i): Fraction(1)} for i in range(n - 1)]

def xx_label(n, i):
    return "".join("X" if k in (i, i + 1) else "I" for k in range(n))

def bounded_closure(controls, seeds, cap=512, work_cap=100000):
    if cap < 1 or work_cap < 0: raise ValueError("cap must be positive and work_cap nonnegative")
    basis = set(seeds)
    if len(basis) > cap: return set(sorted(basis)[:cap]), False, 0, ["INITIAL_CAP"]
    queue = deque(sorted(basis)); operations = 0; omissions = set(); complete = True
    while queue:
        o = queue.popleft()
        for h in controls:
            if operations >= work_cap:
                return basis, False, operations, sorted(omissions | {"WORK_CAP"})
            operations += 1
            for z in commutator_i(h, o):
                if z not in basis:
                    if len(basis) >= cap:
                        omissions.add(z); complete = False
                    else:
                        basis.add(z); queue.append(z)
    return basis, complete, operations, sorted(omissions)

def run(max_n=6, cap=512, work_cap=100000):
    rows = []
    for n in range(2, max_n + 1):
        controls = quadratic_generators(n)
        seeds = [next(iter(h)) for h in controls]
        t = time.perf_counter(); basis, complete, ops, omissions = bounded_closure(controls, seeds, cap, work_cap)
        rows.append({"n": n, "family": "quadratic", "size": len(basis), "expected": n * (2*n - 1),
                     "complete": complete, "operations": ops, "seconds": time.perf_counter() - t,
                     "omissions": omissions[:8], "verified_zero_residual": complete and all(not any(generator_matrix(h, sorted(basis))[2].values()) for h in controls), "cap": cap, "work_cap": work_cap})
        rows[-1]['total_seconds_including_verification'] = time.perf_counter() - t
        controls2 = controls + [{"Z" + "Z" + "I" * (n - 2): Fraction(1)}]
        t = time.perf_counter(); b2, c2, op2, om2 = bounded_closure(controls2, seeds, cap, work_cap)
        rows.append({"n": n, "family": "quadratic_plus_ZZ", "size": len(b2), "complete": c2,
                     "operations": op2, "seconds": time.perf_counter() - t, "omissions": om2[:8],
                     "central_exception": n == 2, "verified_zero_residual": c2 and all(not any(generator_matrix(h, sorted(b2))[2].values()) for h in controls2), "cap": cap, "work_cap": work_cap})
        rows[-1]['total_seconds_including_verification'] = time.perf_counter() - t
    return rows

if __name__ == "__main__":
    out = Path(__file__).parents[1] / "results" / "scaling_results.json"; out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(run(), indent=2) + "\n")
    print(out)
