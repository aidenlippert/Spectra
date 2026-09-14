"""A full completed-spin ceiling proposed by a small-subspace solve."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time

from research.certificate_scaling.commutator_dual_witness import moment_decode, psd, seed
from research.joint_patterns_20260913.dual import affine_round, check as check_baseline, cone_groups, integer_grams
from research.molecular_collective_20260913.core import digest, extract, factor_operators
from research.spin_completion_20260913.core import frames
from research.spin_enrichment_20260913.campaign import ROOT, save


def full_groups(p, tail, masks):
    inherited, signature = cone_groups(p, tail, masks)
    patterns = [q for _, q in factor_operators(p, tail)]
    added = frames(patterns, p['modes'], signature)
    if [len(g['polynomials']) for g in added] != [30, 30, 93, 93, 93, 93, 30, 30]:
        raise ValueError('Full completed-spin frame is incomplete')
    return [g for g in inherited if g['kind'] == 'square'], added


def check(data, tail, witness):
    start = time.monotonic()
    if witness.get('kind') != 'spin_completion_full_dual_v1':
        raise ValueError('Unknown completed-spin full dual')
    # Includes input binding, coefficient box, all ideal and quadratic gates.
    # The inherited 33-dimensional spin-summed PSD gates are redundant inside
    # the complete 93-dimensional frames checked explicitly below.
    inherited = check_baseline(data, tail, {**witness, 'kind': 'joint_density_dual_v1'})
    p = extract(data, tail['center_number']); _, groups = full_groups(p, tail, witness['parity_masks'])
    y = moment_decode(witness['moments'], p['modes']); stats = []
    for group in groups:
        matrices, _ = integer_grams(group['polynomials'], [y], 'anticommutator')
        before = time.monotonic(); result = psd(matrices[0])
        stats.append({**result, 'PSD_seconds': time.monotonic()-before})
    return {'retained_lower_ceiling_Ha': inherited['retained_lower_ceiling_Ha'],
        'original_lower_ceiling_Ha': inherited['original_lower_ceiling_Ha'],
        'original_lower_ceiling_float_Ha': inherited['original_lower_ceiling_float_Ha'],
        'inherited_checks': inherited, 'full_completed_spin_PSD_checks': stats,
        'replay_seconds': time.monotonic()-start, 'many_body_states_enumerated': 0,
        'scope': 'Ceiling on every lower in the entire quadratic plus 492-generator completed-spin frame, all selected subspaces, body-one number ideal, coefficient-L1 residual and fixed lower tail shift.'}


def propose(source, out):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    import numpy as np
    from scipy.linalg import eigh
    prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    raw_path = source/'dual_proposal.json'; raw = json.loads(raw_path.read_text())
    cert = json.loads((source/'certificate.json').read_text()); metadata = json.loads((source.parent/'summary.json').read_text())
    if cert.get('kind') != 'spin_completion_subspace_v1' or cert.get('fixture_sha256') != digest(data) or cert.get('tail_sha256') != digest(tail):
        raise ValueError('Source proposal concerns different frozen inputs')
    masks = metadata['model_construction']['symmetry']['parity_masks']; p = extract(data, tail['center_number'])
    baseline, added = full_groups(p, tail, masks)
    rounded, rank = affine_round(raw, p['modes'], p['particles']); allwords = set(rounded)
    for k in range(3):
        for inds in combinations(range(p['modes']), k):
            allwords.add(tuple((1, i) for i in inds)+tuple((0, i) for i in inds))
    trace = {w: seed(w, p['modes'], p['particles']) for w in allwords}
    thresholds = []
    for group in baseline+added:
        pair, denominator = integer_grams(group['polynomials'], [rounded, trace], group['kind'])
        candidate, anchor = [np.array([[float(F(v, denominator)) for v in row] for row in matrix]) for matrix in pair]
        value = float(eigh(candidate, anchor, eigvals_only=True, subset_by_index=[0, 0])[0])
        threshold = max(0., -value/(1.-value)) if value < 0 else 0.
        thresholds.append({'kind': group['kind'], 'dimension': len(candidate), 'minimum_generalized_moment': value,
            'proposed_mixture_threshold': threshold})
    required = max(r['proposed_mixture_threshold'] for r in thresholds)
    save(out/'proposal.json', {'source': str(source.resolve()), 'affine_rank': rank, 'thresholds': thresholds,
        'proposed_mixture_threshold': required, 'seconds': time.monotonic()-start})
    print(json.dumps({'stage': 'mixture_proposed', 'required': required, 'seconds': time.monotonic()-start}), flush=True)
    attempts = []; tried = set()
    for margin in (1.02, 1.1, 1.5, 2.):
        mix = F(max(10000, math.ceil(required*margin*10**12)), 10**12)
        if mix in tried:
            continue
        tried.add(mix)
        if mix > F(1, 100):
            attempts.append({'mix': str(mix), 'refusal': 'Exceeds 0.01 trace-mixture cap'}); continue
        y = {w: (1-mix)*rounded.get(w, F(0))+mix*trace[w] for w in allwords}
        witness = {'kind': 'spin_completion_full_dual_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
            'parity_masks': masks, 'trace_mixture': str(mix),
            'moments': [{'word': w, 'value': str(v)} for w, v in sorted(y.items(), key=lambda x: (len(x[0]), x[0])) if v or not w]}
        save(out/'pending_witness.json', witness)
        before = time.monotonic()
        try:
            receipt = check(data, tail, witness)
        except ValueError as error:
            attempts.append({'mix': str(mix), 'refusal': str(error), 'seconds': time.monotonic()-before})
            save(out/'attempts.json', attempts); continue
        receipt.update({'attempts': attempts, 'trace_mixture': str(mix), 'source_proposal_sha256': hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            'source_certificate_sha256': hashlib.sha256((source/'certificate.json').read_bytes()).hexdigest(),
            'wall_seconds': time.monotonic()-start, 'affine_rank': rank})
        (out/'witness.json').write_text(json.dumps(witness, separators=(',', ':'))+'\n'); save(out/'receipt.json', receipt)
        print(json.dumps(receipt), flush=True); return receipt
    save(out/'attempts.json', attempts)
    raise ValueError('No exact completed-spin witness accepted within the declared schedule')


def run(source, out):
    log_path = out.with_suffix('.log'); watchdog = out.with_name(out.name+'_watchdog.json')
    if out.exists() or log_path.exists() or watchdog.exists():
        raise FileExistsError('Refusing to overwrite a full-dual diagnostic')
    start = time.monotonic(); failure = None; budget = 300.
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    with log_path.open('w') as log:
        try:
            proc = subprocess.run([sys.executable, '-m', 'research.spin_enrichment_20260913.full_dual',
                '--source', str(source), '--out', str(out), '--worker'], env=environment, cwd=ROOT,
                stdout=log, stderr=subprocess.STDOUT, timeout=budget)
            if proc.returncode:
                failure = f'exit_{proc.returncode}'
        except subprocess.TimeoutExpired:
            failure = 'end_to_end_wall_timeout'
    result = {'budget_seconds': budget, 'wall_seconds': time.monotonic()-start, 'failure': failure,
        'accepted': (out/'receipt.json').exists() and failure is None}
    save(watchdog, result); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.worker:
        propose(args.source, args.out)
    else:
        run(args.source, args.out)
