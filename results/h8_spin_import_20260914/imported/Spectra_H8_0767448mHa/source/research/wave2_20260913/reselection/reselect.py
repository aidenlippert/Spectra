"""Bounded expand, Ritz solve, prune, and reselection experiment."""
from __future__ import annotations
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from research.all_angles_20260913.selected_refinement.refine import CompiledHamiltonian, restricted_matrix, ritz
from research.certificate_scaling.streaming_reference_upper import upper


def _witness(coeffs):
    amps = [int(round(float(x) * 10**12)) for x in coeffs]
    return {'states': [s for s, a in zip(_witness.states, amps) if a],
            'amplitudes': [a for a in amps if a]}


def run(fixture_path, input_path, out_path, sweeps=3, prior_receipt_path=None):
    started = time.monotonic()
    out = Path(out_path); out.mkdir(parents=True, exist_ok=False)
    fixture_raw = Path(fixture_path).read_bytes(); fixture = json.loads(fixture_raw)
    input_raw = Path(input_path).read_bytes(); supplied = json.loads(input_raw)['independent_upper']
    if len(supplied['states']) != 1024 or len(supplied['states']) != len(supplied['amplitudes']):
        raise ValueError('input witness must contain exactly 1024 states and amplitudes')
    upper(fixture, supplied)
    compiled = CompiledHamiltonian(fixture)
    basis = list(supplied['states'])
    import numpy as np
    coeff = np.asarray(supplied['amplitudes'], dtype=float) / 10**12
    coeff /= np.linalg.norm(coeff)
    rows = []
    best = None
    candidate_total = 0

    def record(stage, basis, coeff, numeric, nnz):
        nonlocal best
        _witness.states = basis
        wit = _witness(coeff)
        exact, replay = upper(fixture, wit)
        raw = (json.dumps({'independent_upper': wit, 'upper': str(exact)}, separators=(',', ':'))+'\n').encode()
        name = f'{stage}_upper.json'; (out/name).write_bytes(raw)
        row = {'stage': stage, 'basis_size': len(basis), 'witness_support': len(wit['states']),
               'numeric_ritz': float(numeric), 'upper': str(exact), 'upper_float': float(exact),
               'matrix_nnz': int(nnz),
               'witness_sha256': sha256(raw).hexdigest(), 'witness_file': name,
               'candidate_count_cumulative': candidate_total}
        rows.append(row)
        if best is None or exact < best[0]: best = (exact, name)
        return wit

    # Initial input is retained as the starting point and exact endpoint.
    matrix = restricted_matrix(compiled, basis)
    energy, coeff = ritz(matrix, coeff)
    _witness.states = basis; record('initial_1024', basis, coeff, energy, matrix.nnz)
    for sweep in range(1, sweeps + 1):
        present = set(basis); residual = {}
        for state, c in zip(basis, coeff):
            for target, value in compiled.action(state).items():
                if target not in present: residual[target] = residual.get(target, 0.0) + value*c
        scored = sorted(((abs(v)**2 / max(abs(compiled.diagonal(s)-energy), 1e-3), s)
                         for s, v in residual.items() if v), reverse=True)
        candidate_total += len(residual)
        expanded = sorted(present | {s for _, s in scored[:1024]})
        m = restricted_matrix(compiled, expanded); e, v = ritz(m)
        _witness.states = expanded; record(f'sweep_{sweep:02d}_expanded_2048', expanded, v, e, m.nnz)
        keep = np.argsort(np.abs(v))[-1024:]
        basis = sorted(expanded[i] for i in keep)
        m = restricted_matrix(compiled, basis); energy, coeff = ritz(m)
        _witness.states = basis; record(f'sweep_{sweep:02d}_pruned_1024', basis, coeff, energy, m.nnz)
    prior = json.loads(Path(prior_receipt_path).read_text()) if prior_receipt_path else None
    receipt = {'kind': 'bounded_expand_prune_reselection_v1', 'fixture_sha256': sha256(fixture_raw).hexdigest(),
               'input_sha256': sha256(input_raw).hexdigest(), 'fixture_path': str(fixture_path),
               'input_path': str(input_path), 'sweeps': sweeps, 'max_discovery_support': 2048,
               'final_support': 1024, 'candidate_count_total': candidate_total,
               'hamiltonian_terms': len(compiled.terms), 'full_fock_enumeration': False,
               'best_exact_upper': str(best[0]), 'best_witness_file': best[1], 'rows': rows,
               'target_upper': -9.254455014194086,
               'target_met_by_final_1024': bool(rows[-1]['upper_float'] < -9.254455014194086),
               'prior_run_wall_seconds': prior.get('wall_seconds') if prior else None,
               'prior_run_receipt_sha256': sha256(Path(prior_receipt_path).read_bytes()).hexdigest() if prior_receipt_path else None,
               'total_wall_seconds_including_prior': (prior.get('wall_seconds', 0.0) if prior else 0.0) + time.monotonic()-started,
               'wall_seconds': time.monotonic()-started}
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('--fixture', required=True); p.add_argument('--input', required=True); p.add_argument('--out', required=True); p.add_argument('--prior-receipt'); a=p.parse_args()
    print(json.dumps(run(a.fixture, a.input, a.out, prior_receipt_path=a.prior_receipt), indent=2))
