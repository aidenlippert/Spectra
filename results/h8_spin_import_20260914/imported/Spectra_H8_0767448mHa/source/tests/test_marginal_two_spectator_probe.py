"""Independent bit-action validation of the proposed two-spectator class."""
from pathlib import Path
import importlib.util

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'results/marginal_graded_hubbard8/discovery/two_spectator_overlap_probe.py'
spec=importlib.util.spec_from_file_location('two_spectator_probe',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_every_two_spectator_car_action_matches_independent_bit_swaps():
    for label in module.LABELS:
        data=module.operator(label)
        i,j,k,l,p,q=label
        terms=[(i,j,k,l,p,q,1),(4-j,4-i,4-l,4-k,q,p,-1),
               (i+1,j+1,k+1,l+1,p,q,-1),(5-j,5-i,5-l,5-k,q,p,1)]
        for state in range(4096):
            charges=[((state>>(2*s))&3).bit_count()-1 for s in range(6)]
            row={}
            for a,b,c,d,r,s,weight in terms:
                assert len({a,b,c,d})==4 and r in (1,2) and s in (1,2)
                weight*=charges[c]**r*charges[d]**s
                for spin in (0,1):
                    left,right=2*a+spin,2*b+spin
                    if ((state>>left)&1)==((state>>right)&1):continue
                    sign=(-1)**sum((state>>mode)&1 for mode in range(left+1,right))
                    target=state^(1<<left)^(1<<right)
                    row[target]=row.get(target,0)+weight*sign
            assert data[state]=={s:a for s,a in row.items() if a}


def test_basis_covers_every_parity_allowed_nonconstant_product():
    from itertools import combinations,product
    candidates=set()
    for i,j in combinations(range(5),2):
        for k,l in combinations(set(range(5))-{i,j},2):
            for p,q in product((1,2),repeat=2):
                if (p+q+j-i)%2==1:candidates.add((i,j,k,l,p,q))
    covered=set(module.LABELS)|{module.reflected(label) for label in module.LABELS}
    assert len(candidates)==60 and len(module.LABELS)==30
    assert candidates==covered
    assert module.hopping_norm_bound()==2
