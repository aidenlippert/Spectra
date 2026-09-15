"""Bounded, explicitly enumerative discovery for signed molecular routing."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import time

from research.interference_routing_20260913.routing import (
    load_sector, components, network, positive_tree, replay, check_network, MAX_COMPONENT,
)

ROOT = Path(__file__).resolve().parents[2]


def propose(path, method, amplitude_kind):
    start = time.monotonic()
    import numpy as np
    a, states, data, sha = load_sector(path)
    assembled = time.monotonic()
    proof = {'method': method, 'fixture_sha256': sha, 'modes': data['modes'],
             'particles': data['particles'], 'components': []}
    uppers = []
    diagnostics = []
    eig_work = 0
    for indices in components(a):
        if len(indices) > MAX_COMPONENT:
            raise ValueError('Component work budget')
        block = [[a[i][j] for j in indices] for i in indices]
        ev, u = np.linalg.eigh(np.array(block, dtype=float))
        eig_work += len(block)**3
        v = [int(round(float(x)*10**9)) for x in u[:, 0]]
        uppers.append((float(ev[0]), indices, v))
        phi = v[:] if amplitude_kind == 'eigenvector' else [1]*len(block)
        # Exact zero amplitudes cannot enter a nonsingular congruence.
        if any(x == 0 for x in phi):
            raise ValueError('Rounded trial has nodes; split proposal refused')
        lower = F(math.floor(float(ev[0])*10**8)-1, 10**8) if method == 'grounded_tree_routing_v1' else None
        edges, local = network(block, phi, lower)
        tree = positive_tree(len(block)+(lower is not None), edges)
        item = {'states': [states[i] for i in indices], 'amplitudes': phi, 'tree': tree}
        if lower is not None:
            item['lower'] = str(lower)
        proof['components'].append(item)
        # A finite falsification: try the same amplitude without diagonal slack.
        split_edges, _ = network(block, phi)
        try:
            split_tree = positive_tree(len(block), split_edges)
            check_network(block, phi, split_tree)
            split_status = 'accepted'
        except ValueError as exc:
            split_status = str(exc)
        diagnostics.append({'dimension': len(block), 'numerical_ground': float(ev[0]),
                            'min_local_energy': str(min(local)), 'split_status': split_status})
    _, indices, v = min(uppers, key=lambda x: x[0])
    proof['upper_amplitudes'] = [0]*len(a)
    for i, x in zip(indices, v):
        proof['upper_amplitudes'][i] = x
    generated = time.monotonic()
    receipt = replay(path, proof)
    proof['claim'] = {k: receipt[k] for k in ('lower', 'upper', 'width')}
    accounting = {'assembly_seconds': assembled-start, 'proposal_and_diagnostics_seconds': generated-assembled,
                  'total_seconds_including_replay': time.monotonic()-start,
                  'dense_eigensolve_cubic_work_proxy': eig_work,
                  'amplitude_kind': amplitude_kind, 'diagnostics': diagnostics,
                  'discovery': 'Full fixed-N CAR matrix and dense eigensolves in all connected components'}
    return proof, receipt, accounting


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    cases = [
        ('h4_line_guided', 'results/certificate_scaling/active_space_ladder/h4/fixture.json', 'split_tree_routing_v1', 'eigenvector'),
        ('h4_line_grounded', 'results/certificate_scaling/active_space_ladder/h4/fixture.json', 'grounded_tree_routing_v1', 'uniform'),
        ('h4_rectangle_grounded', 'results/marginal_molecule/h4_rectangle_sto3g.json', 'grounded_tree_routing_v1', 'uniform'),
        ('h4_square_grounded', 'results/marginal_molecule_stress/h4_square_sto3g.json', 'grounded_tree_routing_v1', 'uniform'),
        ('h6_budget', 'results/certificate_scaling/active_space_ladder/h6/fixture.json', 'grounded_tree_routing_v1', 'uniform'),
    ]
    results = []
    for name, relative, method, kind in cases:
        path = ROOT/relative
        started = time.monotonic()
        try:
            proof, receipt, accounting = propose(path, method, kind)
            (out/(name+'_proof.json')).write_text(json.dumps(proof, separators=(',', ':'))+'\n')
            result = {'name': name, 'fixture': str(path), 'status': 'accepted',
                      'proof_bytes': (out/(name+'_proof.json')).stat().st_size,
                      'receipt': receipt, 'accounting': accounting}
        except ValueError as exc:
            result = {'name': name, 'fixture': str(path), 'status': 'refused', 'reason': str(exc),
                      'seconds': time.monotonic()-started}
        results.append(result)
        print(json.dumps({'name': name, 'status': result['status'], 'width': result.get('receipt', {}).get('width_float'),
                          'reason': result.get('reason')}), flush=True)
    # Recompute the actual failed baseline with the unchanged accepting checker.
    from experiments.marginal_symbolic import verified_residual, decode
    source = ROOT/'results/certificate_scaling/commutator_dictionary/h4_creator_channels/certificate.json'
    raw = source.read_bytes()
    cert = json.loads(raw)
    fixture = json.loads((ROOT/cases[0][1]).read_text())
    if cert['modes'] != fixture['modes'] or cert['particles'] != fixture['particles'] or decode(cert['hamiltonian'], 8, 4) != decode(fixture['hamiltonian'], 8, 4):
        raise ValueError('Baseline is a different physical problem')
    baseline_start = time.monotonic()
    _, old = verified_residual(cert)
    first = results[0]['receipt']
    baseline = {'certificate': str(source), 'sha256': hashlib.sha256(raw).hexdigest(),
                'lower': old['lower'], 'upper': first['upper'],
                'width': str(F(first['upper'])-F(old['lower'])),
                'width_float': float(F(first['upper'])-F(old['lower'])),
                'replay_seconds': time.monotonic()-baseline_start,
                'comparison': 'Same rational H and exact upper; routing changes the proof family and enumerates every state.'}
    payload = {'baseline': baseline, 'cases': results}
    (out/'summary.json').write_text(json.dumps(payload, indent=2)+'\n')
    return payload


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    run(p.parse_args().out)
