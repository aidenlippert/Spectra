#!/bin/sh
set -eu
ssh -o BatchMode=yes -o ServerAliveInterval=20 -i /Users/aidenlippert/.ssh/id_ed25519 ubuntu@129.146.98.55 'cd /home/ubuntu/spectra-gpu-20260912 && CUDA_PATH=/usr OPENBLAS_NUM_THREADS=1 timeout 300 /home/ubuntu/spectra-gpu-venv/bin/python resident_batch_benchmark.py --data resident.npz --out resident_recheck.json --evaluations 1024 --batch 256'
