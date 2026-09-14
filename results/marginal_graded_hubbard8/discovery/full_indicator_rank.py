"""Exact rank of the full binary-indicator basis and actual correction span."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_quadratic_charge_telescope import LABELS as Q,local_value as qvalue
from experiments.marginal_charge_square_pairs import LABELS as S,local_value as svalue
from experiments.marginal_charge_indicator_telescope import LABELS as I,ALL_LABELS,five_site_value,local_value as ivalue


def rank(rows):
    pivots={}
    for row in rows:
        row=list(map(F,row))
        for j in sorted(pivots):
            a=row[j]
            if a:row=[v-a*b for v,b in zip(row,pivots[j])]
        j=next((j for j,v in enumerate(row) if v),None)
        if j is not None:
            a=row[j];pivots[j]=[v/a for v in row]
    return len(pivots)


def main():
    source=ROOT/'results/marginal_graded_hubbard8/charge_square_pairs/W_zero/range_two_family_limit_certificate.json'
    c=json.loads(source.read_text())
    shapes=[{int(s):F(v) for s,v in d.items()} for d in c['diagonal_shapes']]
    indicator_rank=rank([five_site_value(sum((0 if mask>>i&1 else 1)*4**i for i in range(5)),{key:F(1)}) for key in ALL_LABELS] for mask in range(32))
    combined_rank=rank([d.get(s&1023,F(0))-d.get(s>>2,F(0)) for d in shapes]+[fn(s,{key:F(1)}) for names,fn in [(Q,qvalue),(S,svalue),(I,ivalue)] for key in names] for s in range(4096))
    if indicator_rank!=12 or combined_rank!=25:raise ValueError('Unexpected exact indicator or correction rank')
    counts=Counter(map(len,ALL_LABELS.values()))
    if counts!={1:2,2:4,3:4,4:2}:raise ValueError('Incomplete reflection-odd Boolean monomial basis')
    files={Path(__file__).resolve(),source}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'binary_patterns':32,'indicator_rank':indicator_rank,'indicator_basis_counts':dict(counts),'full_fock_determinants':4096,'actual_diagonal_span_rank':combined_rank,'actual_diagonal_columns':25,'higher_indicator_labels':list(I),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Exact finite rank on the full five-site binary space and six-site Fock space. The additional six indicator directions are independent of the actual nine sparse shapes, six quadratic and four square-pair directions. Does not prove global representability.'}
    (ROOT/'results/marginal_graded_hubbard8/full_charge_indicators/full_indicator_rank.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__=='__main__':main()
