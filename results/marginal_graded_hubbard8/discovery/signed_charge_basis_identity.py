"""Exact equivalence of the disjoint charge patterns and full polynomial space."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
from math import prod
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_signed_charge_telescope import PATTERNS,ORBIT,five_site_value


def main():
    charges=list(product((-1,0,1),repeat=5))
    fixed_ph=sum(q==tuple(-v for v in q) for q in charges)
    fixed_reflection=sum(q==q[::-1] for q in charges)
    fixed_both=sum(q==tuple(-v for v in q[::-1]) for q in charges)
    dimension=(243+fixed_ph-fixed_reflection-fixed_both)//4
    if dimension!=52 or len(PATTERNS)!=52 or len(ORBIT)!=208:raise ValueError('Incomplete symmetry-projected basis')
    exponents=[e for e in product(range(3),repeat=5) if sum(e)%2==0 and e<e[::-1]]
    if len(exponents)!=52:raise ValueError('Incomplete polynomial basis')
    def value(q,e):return prod(v**a for v,a in zip(q,e))-prod(v**a for v,a in zip(q,e[::-1]))
    for e in exponents:
        terms={key:F(value(q,e)) for key,q in PATTERNS.items()}
        for s in range(1024):
            q=tuple(((s>>(2*i))&3).bit_count()-1 for i in range(5))
            if five_site_value(s,terms)!=value(q,e):raise ValueError('Polynomial and pattern operator differ')
    files={Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'signed_charge_patterns':243,'fixed_ph_patterns':fixed_ph,'fixed_reflection_patterns':fixed_reflection,'fixed_ph_reflection_patterns':fixed_both,'projected_dimension':dimension,'disjoint_nonzero_patterns':208,'polynomials_reconstructed':52,'physical_determinants_per_polynomial':1024,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Exact equivalence on physical states of the complete PH-even reflection-odd charge-polynomial and disjoint pattern bases. No numerical conditioning improvement or quantum representability is inferred from this identity.'}
    (ROOT/'results/marginal_graded_hubbard8/full_signed_charge/basis_identity.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__=='__main__':main()
