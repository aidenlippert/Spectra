"""Exact missing diagonal-consistency moments in the accepted local dual."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_projector_extendibility import _telescoping_diagonal


def main():
    started = time.monotonic()
    base = ROOT/'results/marginal_graded_hubbard8/joint_projector/signed_density'
    folder = base/'symmetric_diagonals_9/extra_31'
    certificate_path = folder/'diagonal_family_limit_certificate.json'
    receipt_path = folder/'diagonal_family_limit_replay.json'
    shapes_path = base/'symmetry_diagonal_candidates.json'
    receipt = json.loads(receipt_path.read_text())
    if not receipt['accepted'] or receipt['source_sha256'][str(certificate_path.relative_to(ROOT))] != hashlib.sha256(certificate_path.read_bytes()).hexdigest():
        raise ValueError('Mixture must match the accepted exact receipt')
    c = json.loads(certificate_path.read_text())
    mismatch = [F(0)]*1024
    for state in c['mixture']:
        vector = {int(s): v for s, v in state['vector'].items()}
        norm = sum(v*v for v in vector.values())
        weight = F(state['weight'])
        for s, v in vector.items():
            p = weight*F(v*v, norm)
            mismatch[s & 1023] += p
            mismatch[s >> 2] -= p
    candidates = json.loads(shapes_path.read_text())['candidates']
    moments = []
    for i, item in enumerate(candidates):
        shape = _telescoping_diagonal(item['diagonal'])
        moment = sum((v*mismatch[s] for s, v in shape.items()), F(0))
        moments.append({'basis_index': i, 'exact_moment': str(moment),
                        'moment_float': float(moment), 'entries': len(shape)})
    selected = list(range(8))+[31]
    if any(F(moments[i]['exact_moment']) for i in selected):
        raise ValueError('A supposedly canceled family moment is nonzero')
    violations = sorted((row for row in moments if F(row['exact_moment'])), key=lambda row: abs(F(row['exact_moment'])), reverse=True)
    files = [Path(__file__).resolve(), certificate_path, receipt_path, shapes_path,
             ROOT/'experiments/marginal_projector_extendibility.py']
    result = {'accepted': True, 'tested_shapes': len(moments), 'selected_shapes_cancel_exactly': selected,
              'violated_shapes': len(violations), 'largest_violations': violations[:10],
              'moments': moments, 'seconds': time.monotonic()-started,
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'scope': 'Exact necessary translation-consistency violations in the accepted local PSD dual mixture. Any nonzero moment of Y_left-Y_right prevents it from being a translation-invariant physical marginal. These are separators, not stronger accepted energy bounds. The current 64-entry cap prevents simply adding every remaining shape.'}
    (folder/'remaining_diagonal_overlap_moments.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('moments', 'source_sha256', 'largest_violations')}, indent=2))
    print(json.dumps([{k: v for k, v in row.items() if k != 'exact_moment'} for row in violations[:10]]), flush=True)


if __name__ == '__main__':
    main()
