"""Exact ceiling for the complete spin-resolved frame, not just selected cuts."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

from research.certificate_scaling.commutator_dual_witness import moment_decode, psd, seed
from research.joint_patterns_20260913.dual import affine_round, check as check_spin_summed, cone_groups, integer_grams
from research.joint_patterns_20260913.spin_diagnostic import spin_generator
from research.molecular_collective_20260913.core import digest, extract, factor_operators


def frames(p, tail, masks):
    old_groups, signature = cone_groups(p, tail, masks)
    patterns = [q for _, q in factor_operators(p, tail)]; parts = {}
    for k, spin in [(-1, -1)]+[(k, spin) for k in range(10) for spin in range(2)]:
        for mode in range(p['modes']):
            q = spin_generator(patterns, [k, spin, mode], p['modes'])
            keys = {signature(w) for w in q}
            if not q or len(keys) != 1:
                raise ValueError('Full spin generator lacks definite symmetry')
            parts.setdefault(next(iter(keys)), []).append(q)
    added = [{'kind': 'anticommutator', 'polynomials': parts[key]} for key in sorted(parts)]
    if sum(len(g['polynomials']) for g in added) != 252:
        raise AssertionError('Full spin frame is incomplete')
    return [g for g in old_groups if g['kind'] == 'square'], added


def check(data, tail, witness):
    start = time.monotonic()
    if witness.get('kind') != 'full_spin_frame_dual_v1':
        raise ValueError('Unknown full-spin witness')
    # This reuses all binding, normalization, Hermiticity, coefficient-box,
    # symmetry, ideal, and baseline PSD gates. Its four old density blocks
    # are redundant positive subspaces of the full spin frame checked below.
    base = check_spin_summed(data, tail, {**witness, 'kind': 'joint_density_dual_v1'})
    p = extract(data, tail['center_number']); _, added = frames(p, tail, witness['parity_masks'])
    y = moment_decode(witness['moments'], p['modes']); stats = []
    for group in added:
        matrices, _ = integer_grams(group['polynomials'], [y], 'anticommutator')
        stats.append(psd(matrices[0]))
    return {'retained_lower_ceiling_Ha': base['retained_lower_ceiling_Ha'],
        'original_lower_ceiling_Ha': base['original_lower_ceiling_Ha'],
        'original_lower_ceiling_float_Ha': base['original_lower_ceiling_float_Ha'],
        'inherited_checks': base, 'full_spin_PSD_checks': stats,
        'replay_seconds': time.monotonic()-start, 'many_body_states_enumerated': 0,
        'scope': 'Ceiling for the complete quadratic plus spin-resolved pattern anticommutator frame and all its subspaces, with the body-one number ideal, coefficient-L1 residual, and fixed lower tail shift.'}


def propose(data, tail, source, out):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    raw_path = source/'dual_proposal.json'; raw = json.loads(raw_path.read_text())
    construction = json.loads((source/'construction.json').read_text())
    if construction['directions'] != 252 or construction['separate'] or construction['anti_dimensions'] != [63]*4:
        raise ValueError('Source was not the complete joint spin-frame control')
    cert = json.loads((source/'certificate.json').read_text())
    if cert['fixture_sha256'] != digest(data) or cert['tail_sha256'] != digest(tail):
        raise ValueError('Full-frame source is for different frozen inputs')
    policy = json.loads((source.parent/'summary.json').read_text())['model_construction']
    masks = policy['symmetry']['parity_masks']; p = extract(data, tail['center_number'])
    baseline, added = frames(p, tail, masks)
    rounded, rank = affine_round(raw, p['modes'], p['particles']); allwords = set(rounded)
    for k in range(3):
        for inds in combinations(range(p['modes']), k):
            allwords.add(tuple((1, i) for i in inds)+tuple((0, i) for i in inds))
    trace = {w: seed(w, p['modes'], p['particles']) for w in allwords}
    matrices = [integer_grams(g['polynomials'], [rounded, trace], g['kind'])[0] for g in baseline+added]
    attempts = []
    for mix in (F(0), F(1, 10**8), F(1, 10**7), F(1, 10**6), F(1, 10**5), F(1, 10**4), F(1, 1000)):
        try:
            for candidate, seed_matrix in matrices:
                psd([[(mix.denominator-mix.numerator)*a+mix.numerator*b for a, b in zip(row, sr)] for row, sr in zip(candidate, seed_matrix)])
        except ValueError as error:
            attempts.append({'mix': str(mix), 'refusal': str(error)}); continue
        y = {w: (1-mix)*rounded.get(w, F(0))+mix*trace[w] for w in allwords}
        witness = {'kind': 'full_spin_frame_dual_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
            'parity_masks': masks, 'trace_mixture': str(mix),
            'moments': [{'word': w, 'value': str(v)} for w, v in sorted(y.items(), key=lambda item: (len(item[0]), item[0])) if v or not w]}
        receipt = check(data, tail, witness)
        receipt.update({'attempts': attempts, 'trace_mixture': str(mix), 'affine_rank': rank,
            'source_proposal_sha256': hashlib.sha256(raw_path.read_bytes()).hexdigest(), 'wall_seconds': time.monotonic()-start})
        (out/'witness.json').write_text(json.dumps(witness, separators=(',', ':'))+'\n')
        (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n'); print(json.dumps(receipt), flush=True)
        return receipt
    (out/'refusals.json').write_text(json.dumps(attempts, indent=2)+'\n')
    raise ValueError('No accepted exact full-spin dual within the mixture schedule')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--fixture', type=Path, required=True); parser.add_argument('--tail', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); args = parser.parse_args()
    propose(json.loads(args.fixture.read_text()), json.loads(args.tail.read_text()), args.source, args.out)
