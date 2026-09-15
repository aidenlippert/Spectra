"""Execute full-number-X control by changing only the multiplier Sz predicate."""
from pathlib import Path
import json
from experiments.marginal_symbolic import decode
import research.certificate_scaling.spin_invariant_discovery as backend

ROOT=Path(__file__).resolve().parents[3]
fixture=ROOT/'results/certificate_scaling/spin_irrep/h4/symmetric_hamiltonian.json'
out=ROOT/'results/all_angles_20260913/spin_ablation/h4_full_x'

# Narrow adapter: decomposition, parity, Gram maps, scaling and verifier are
# untouched. Only allowed(w)'s SU(2) predicate is neutralized; parity remains.
backend.spin_weight2=lambda w: 0
f=json.loads(fixture.read_text())
h=decode(f['hamiltonian'],f['modes'],4)
result=backend.run(h,f['modes'],f['particles'],out,solver='SCS',seconds=25)
print(json.dumps(result,indent=2,default=str))
