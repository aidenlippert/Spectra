"""Matched SU(2) multiplier ablation.

The existing invariant discovery certificate is replayed with (a) its twirled
number ideal and (b) the full number conserving degree-2 ideal counted as a
strict superset.  This deliberately separates cone representation from a
solver claim: no projected Hamiltonian or sector-restricted verifier is used.
"""
from pathlib import Path
import json, sys
from experiments.marginal_symbolic import multiplier_basis, decode
from research.certificate_scaling.spin_basis import spin_weight2
from research.certificate_scaling.spin_parity import spin_parity

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results/certificate_scaling/spin_irrep"

def run():
    rows=[]
    for name in ("h4_discovery", "h6_scs"):
        c=json.loads((SRC/name/"certificate.json").read_text())
        m,n=c["modes"],c["particles"]
        h=decode(c["hamiltonian"],m,4); masks=spin_parity(h,m)
        all_basis=[p for p in multiplier_basis(m,2)
                   if all(sum((mask>>i)&1 for _,i in w)%2==0 for mask in masks for w in p)]
        inv=[p for p in all_basis if all(spin_weight2(w)==0 for w in p)]
        twirl_terms=len(decode(c["number_multiplier"],m,4))
        rows.append({"fixture":name,"modes":m,"particles":n,
                     "full_number_conserving_basis":len(all_basis),
                     "restricted_Sz_X_basis":len(inv),
                     "full_number_X_basis":len(all_basis),
                     "larger_basis_terms":len(all_basis)-len(inv),
                     "source_certificate_multiplier_terms":twirl_terms,
                     "original_H_replay":True,
                     "sector_restriction":False,
                     "interpretation":"Sz-weight-zero X versus full number-conserving X; SOS blocks remain SU(2)-invariant. L1 objective is not SU(2)-invariant."})
    out=Path(__file__).resolve().parents[3]/"results/all_angles_20260913/spin_ablation"
    out.mkdir(parents=True,exist_ok=True)
    (out/"receipt.json").write_text(json.dumps({"experiment":"matched_spin_multiplier_ablation","rows":rows},indent=2)+"\n")
    print(json.dumps(rows,indent=2))
if __name__=="__main__": run()
