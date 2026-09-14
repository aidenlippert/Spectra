"""Exact separators of the previous full-spin dual; no energy inference."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time

from research.certificate_scaling.commutator_dual_witness import evaluate, moment_decode
from research.joint_patterns_20260913.core import anticommutator
from research.joint_patterns_20260913.dual import cone_groups, integer_grams
from research.molecular_collective_20260913.core import digest, extract, factor_operators
from research.spin_subspace_20260913.core import combine
from research.spin_subspace_20260913.full_dual import check as check_prior
from research.spin_completion_20260913.core import frames, generator


def check_separator(data, tail, dual, cert):
    start = time.monotonic()
    if cert.get('kind') != 'spin_completion_separator_v1' or cert.get('fixture_sha256') != digest(data) or cert.get('tail_sha256') != digest(tail) or cert.get('dual_sha256') != digest(dual):
        raise ValueError('Completed-spin separator binding failed')
    p = extract(data, tail['center_number']); patterns = [q for _, q in factor_operators(p, tail)]
    refs = cert['generators']; vector = cert['vector']
    if not refs or len(refs) > p['modes']*(4*len(patterns)+1) or len(vector) != len(refs) or any(type(c) is not int for c in vector) or not any(vector):
        raise ValueError('Invalid bounded integer separator')
    q = combine([generator(patterns, ref, p['modes']) for ref in refs], vector, 1)
    if not q:
        raise ValueError('Zero separator operator')
    poly = anticommutator(q, q); y = moment_decode(dual['moments'], p['modes']); value = evaluate(poly, y)
    if value >= 0:
        raise ValueError('This direction does not separate the previous full-spin witness')
    norm = sum(c*c for c in vector)
    return {'exact_expectation': str(value), 'coefficient_norm_squared': str(norm),
        'normalized_expectation': str(value/norm), 'normalized_expectation_float': float(value/norm),
        'operator_monomials': len(q), 'anticommutator_monomials': len(poly),
        'max_degree': max(map(len, poly), default=0), 'replay_seconds': time.monotonic()-start,
        'scope': 'Necessary-condition violation only. Energy gain requires an independent bound.'}


def run(data, tail, dual, out):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    import numpy as np
    previous = check_prior(data, tail, dual)
    print(json.dumps({'stage': 'prior_dual_accepted', 'seconds': time.monotonic()-start}), flush=True)
    p = extract(data, tail['center_number']); patterns = [q for _, q in factor_operators(p, tail)]
    _, signature = cone_groups(p, tail, dual['parity_masks']); groups = frames(patterns, p['modes'], signature)
    y = moment_decode(dual['moments'], p['modes']); proposals = []; stats = []
    for gid, group in enumerate(groups):
        before = time.monotonic(); matrices, denominator = integer_grams(group['polynomials'], [y], 'anticommutator')
        gram = np.array([[float(F(c, denominator)) for c in row] for row in matrices[0]])
        values, vectors = np.linalg.eigh(gram); words = sorted({w for q in group['polynomials'] for w in q})
        C = np.array([[float(q.get(w, 0)) for q in group['polynomials']] for w in words])
        rank = int(np.linalg.matrix_rank(C, tol=1e-10))
        accepted = []
        for j, value in enumerate(values[:2]):
            if value >= -1e-6:
                continue
            vector = list(map(int, np.rint(vectors[:, j]*10**10).astype(np.int64)))
            cert = {'kind': 'spin_completion_separator_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
                'dual_sha256': digest(dual), 'generators': group['generators'], 'vector': vector}
            receipt = check_separator(data, tail, dual, cert)
            path = out/f'separator_{gid}_{j}.json'; path.write_text(json.dumps(cert, separators=(',', ':'))+'\n')
            accepted.append({'group': gid, 'vector': vector, 'proposed_moment': float(value), 'separator': path.name, 'exact': receipt})
        row = {'group': gid, 'name': group['name'], 'dimension': len(values), 'numerical_operator_rank': rank,
            'minimum_proposed_moment': float(values[0]), 'accepted_separators': len(accepted), 'seconds': time.monotonic()-before}
        stats.append(row); proposals.append(accepted); print(json.dumps(row), flush=True)
    result = {'previous_dual': previous, 'frames': stats, 'candidates': proposals,
        'fixture_sha256': digest(data), 'tail_sha256': digest(tail), 'dual_sha256': digest(dual),
        'wall_seconds': time.monotonic()-start, 'many_body_states_enumerated': 0}
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n'); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); root = Path(__file__).resolve().parents[2]
    prior = root/'results/molecular_collective_20260913/campaign/h6'
    run(json.loads((prior/'fixture.json').read_text()), json.loads((prior/'rank_10/tail.json').read_text()),
        json.loads((root/'results/spin_subspace_20260913/full_dual/witness.json').read_text()), args.out)
