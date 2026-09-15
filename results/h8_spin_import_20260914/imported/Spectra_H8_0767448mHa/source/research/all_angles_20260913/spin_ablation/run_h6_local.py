"""Fallback H6 full-number-X run when HOST B is unavailable."""
from pathlib import Path
import json
from experiments.marginal_symbolic import decode
import research.certificate_scaling.spin_invariant_discovery as backend
ROOT=Path(__file__).resolve().parents[3]
f=json.loads((ROOT/'results/certificate_scaling/spin_irrep/h6/symmetric_hamiltonian.json').read_text())
backend.spin_weight2=lambda w: 0
backend.run(decode(f['hamiltonian'],f['modes'],4),f['modes'],f['particles'],ROOT/'results/all_angles_20260913/spin_ablation/h6_full_number_X',solver='SCS',seconds=120)
