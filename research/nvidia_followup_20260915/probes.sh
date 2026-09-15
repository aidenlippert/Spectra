set -eu
cd /home/ubuntu/spectra-nvidia
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/home/ubuntu/.local/bin/uv pip install --python .venv/bin/python sparseqr
/home/ubuntu/.local/bin/uv pip freeze --python .venv/bin/python > installed.txt
.venv/bin/python -B -m research.nvidia_followup_20260915.probe_jobs
