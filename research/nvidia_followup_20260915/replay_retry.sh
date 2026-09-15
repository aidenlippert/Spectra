set -eu
cd /home/ubuntu/spectra-nvidia
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/usr/bin/time -v -o results/nvidia_followup_20260915/flint_transition_test_time.txt .venv/bin/python -B -m unittest research.nvidia_followup_20260915.test_flint_transition research.nvidia_followup_20260915.test_sector_gate
/usr/bin/time -v -o results/nvidia_followup_20260915/h10_compiled_retry_time.txt .venv/bin/python -B -m research.nvidia_followup_20260915.strict_replay results/transfer_solver_20260915/adaptive/h10_correlated_guide sparse_cudss_adaptive results/transfer_solver_20260915/adaptive/h10_correlated_guide/sparse_cudss_adaptive/exact_compiled_retry --rotated results/transfer_solver_20260915/rotated_h10 --compiled-rationals
