"""Finite-response size diagnostic; does not certify a ground energy."""
from experiments.marginal_defect_dicke import DefectDicke
from experiments.marginal_symbolic import add, mono, adj
from pathlib import Path
import json
import time

out = Path(__file__).with_name('local_action_scaling.json')
if out.exists():
    raise ValueError('Preserve previous diagnostic')
rows = []
for m in (5, 10, 20, 50, 100):
    start = time.monotonic()
    model = DefectDicke(m)
    z = model.atom(rights=m // 2)
    p = mono(((1, 0), (0, m + 1)))
    w = model.apply_polynomial(add(p, adj(p)), z)
    norms, max_atoms = [], 0
    for k in range(8):
        norms.append(str(model.inner(w, w)))
        max_atoms = max(max_atoms, len(w))
        w = model.reference_action(w)
    rows.append({'pairs': m, 'modes': 2 * m, 'initial_dicke_rights': m // 2,
                 'h0_applications': 7, 'max_atoms': max_atoms,
                 'max_active_pairs': max(len(a[0]) for a in w), 'squared_norms': norms,
                 'elapsed_seconds': time.monotonic() - start,
                 'scope': 'Single fixed two-pair Hermitian hopping action and seven H0 powers; no energy certificate'})
    print({k: v for k, v in rows[-1].items() if k != 'squared_norms'}, flush=True)
out.write_text(json.dumps(rows, indent=2) + '\n')
