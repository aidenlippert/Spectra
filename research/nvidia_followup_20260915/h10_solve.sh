set -eu
cd /home/ubuntu/spectra-nvidia
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p results/nvidia_followup_20260915
/usr/bin/time -v -o results/nvidia_followup_20260915/h10_hybrid_time.txt .venv/bin/python -B -m research.gpu_acceleration_20260915.solve results/transfer_solver_20260915/adaptive/h10_correlated_guide hybrid_adaptive --seconds 1100 --mu 2 --backend hybrid --adaptive
PYTHONPATH=. .venv/bin/python -B -S - <<'PY'
import json,hashlib
from pathlib import Path
root=Path('results/transfer_solver_20260915/adaptive/h10_correlated_guide')
source=root/'hybrid_adaptive/export/certificate.json'
(root/'stage2/export').mkdir(parents=True,exist_ok=False)
(root/'stage2/export/certificate.json').write_bytes(source.read_bytes())
(root/'selected_proposal.json').write_text(json.dumps({'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'no_new_optimization':True},indent=2)+'\n')
PY
/usr/bin/time -v -o results/nvidia_followup_20260915/h10_replay_time.txt .venv/bin/python -B -S -m research.transfer_solver_20260915.actions results/transfer_solver_20260915/adaptive/h10_correlated_guide replay
