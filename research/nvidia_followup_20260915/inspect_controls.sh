set -eu
cd /home/ubuntu/spectra-nvidia
.venv/bin/python - <<'PY'
import json
from pathlib import Path
for label,path in [('H8',Path('results/transfer_solver_20260915/cases/h8_cold')),('H10',Path('results/transfer_solver_20260915/adaptive/h10_correlated_guide'))]:
    rows=[]
    for tag in ('dense_200','sparse_cpu_200','sparse_hybrid_200','sparse_cudss_200'):
        d=json.loads((path/tag/'discovery.json').read_text())
        h=json.loads((path/tag/'history.json').read_text())[-1]
        rows.append({'tag':tag,'seconds':d['seconds'],'setup':d['setup_seconds'],
            'QR':d['ideal_projection']['seconds'],'rank':d['ideal_projection']['numerical_ideal_rank'],
            'kernel':d['kernel_seconds'],'iterations':d['iterations_completed'],
            'b':h['b'],'unverified_width_mHa':h['unverified_width_mHa']})
    print(label,json.dumps(rows))
PY
tail -3 results/nvidia_followup_20260915/h10_sparse_cudss_adaptive.log
