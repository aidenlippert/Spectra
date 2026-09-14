"""Exact independence of compact quadratic and charge-square pair telescopes."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_quadratic_charge_telescope import LABELS as QUADRATIC,local_value as quadratic
from experiments.marginal_charge_square_pairs import LABELS as SQUARES,local_value as square,five_site_value


def main():
    pivots={}
    for s in range(4096):
        row=[quadratic(s,{label:F(1)}) for label in QUADRATIC]+[square(s,{label:F(1)}) for label in SQUARES]
        for j in sorted(pivots):
            a=row[j];row=[v-a*b for v,b in zip(row,pivots[j])]
        j=next((j for j,v in enumerate(row) if v),None)
        if j is not None:
            a=row[j];pivots[j]=[v/a for v in row]
    if len(pivots)!=10:raise ValueError('New charge-square pair directions are not independent of the quadratic span')
    rows=[]
    for label in SQUARES:
        rows.append({'pair':label,'five_site_nonzero_entries':sum(bool(five_site_value(s,{label:F(1)})) for s in range(1024)),
                     'local_operator_norm':int(max(abs(square(s,{label:F(1)})) for s in range(4096)))})
    files={Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
    result={'accepted':True,'full_fock_determinants':4096,'quadratic_rank':6,'combined_rank':10,'new_independent_directions':4,'pairs':rows,
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
      'scope':'Exact finite operator rank: four charge-square pair telescopes are independent of the complete quadratic charge telescope span. Does not claim global representability.'}
    (ROOT/'results/marginal_graded_hubbard8/charge_square_pairs/charge_square_pair_rank.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__=='__main__':main()
