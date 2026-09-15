"""Independent bit-action validation of three-spectator CAR consistency probes."""
from pathlib import Path
import importlib.util
from itertools import combinations, product

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / 'results/marginal_graded_hubbard8/discovery/three_spectator_overlap_probe.py'
spec = importlib.util.spec_from_file_location('three_spectator_probe', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_all_three_spectator_car_actions_match_independent_bit_swaps():
    for label in module.LABELS:
        data = module.operator(label)
        i, j, k, l, m, p, q, r = label
        terms = [(i,j,k,l,m,p,q,r,1), (4-j,4-i,4-m,4-l,4-k,r,q,p,-1),
                 (i+1,j+1,k+1,l+1,m+1,p,q,r,-1), (5-j,5-i,5-m,5-l,5-k,r,q,p,1)]
        for state in range(4096):
            charges = [((state >> (2*s)) & 3).bit_count() - 1 for s in range(6)]
            row = {}
            for a,b,c,d,e,u,v,w,weight in terms:
                assert len({a,b,c,d,e}) == 5
                weight *= charges[c]**u * charges[d]**v * charges[e]**w
                for spin in (0,1):
                    left, right = 2*a+spin, 2*b+spin
                    if ((state >> left) & 1) == ((state >> right) & 1):
                        continue
                    sign = (-1)**sum((state >> mode) & 1 for mode in range(left+1,right))
                    target = state ^ (1 << left) ^ (1 << right)
                    row[target] = row.get(target,0) + weight*sign
            assert data[state] == {s:a for s,a in row.items() if a}


def test_three_spectator_basis_covers_all_nonzero_reflection_pairs():
    candidates = set()
    for i,j in combinations(range(5),2):
        k,l,m = sorted(set(range(5)) - {i,j})
        for p,q,r in product((1,2), repeat=3):
            if (p+q+r+j-i) % 2 == 1:
                candidates.add((i,j,k,l,m,p,q,r))
    self_reflected = {x for x in candidates if module.reflected(x) == x}
    covered = set(module.LABELS) | {module.reflected(x) for x in module.LABELS}
    assert len(candidates) == 40 and len(self_reflected) == 4
    assert len(module.LABELS) == 18 and covered == candidates - self_reflected
    for label in self_reflected:
        assert all(not row for row in module.operator(label))
    assert module.hopping_norm_bound() == 2
