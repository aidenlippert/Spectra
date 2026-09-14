"""Propose untrusted speed witnesses, or independently replay them with integers.

Proposal mode intercepts PSD calls solely to collect matrices and propose R;
its output is explicitly nonaccepting. Replay mode never installs that hook,
requires only the standard library, and reconstructs the entire certificate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_congruence_psd import matrix_digest


def propose(certificate, output):
    import numpy as np
    from scipy.linalg import solve_triangular
    import experiments.marginal_projector_extendibility as module
    started = time.monotonic()
    witnesses = {}
    dimensions = []
    failures = []

    def collect(matrix, initial_divisor=1):
        n = len(matrix)
        dimensions.append(n)
        key = matrix_digest(matrix)
        if key not in witnesses and n >= 20:
            try:
                scale = max(matrix[i][i] for i in range(n))
                a = np.array(matrix, dtype=float)/float(scale)
                factor = solve_triangular(np.linalg.cholesky(a).T, np.eye(n), lower=False)
                witnesses[key] = [[round(float(v)*10**12) for v in row] for row in factor]
            except (np.linalg.LinAlgError, ValueError, OverflowError, ZeroDivisionError) as error:
                failures.append({'matrix_sha256': key, 'dimension': n, 'error': str(error)})
        # Shape metadata only, intentionally NOT a PSD acceptance receipt.
        return {'dimension': n, 'untrusted_proposal_collection': True}

    original = module.integer_psd
    try:
        module.integer_psd = collect
        replay(certificate)
    finally:
        module.integer_psd = original
    data = {'accepted': False, 'kind': 'untrusted_integer_congruence_proposal_v1',
            'witnesses': witnesses, 'blocks_collected': len(dimensions),
            'dimensions': dimensions, 'proposal_failures': failures,
            'seconds': time.monotonic()-started,
            'scope': 'Numerical change-of-basis proposals only. Every PSD claim requires fresh exact replay.'}
    output.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    print(json.dumps({k: v for k, v in data.items() if k not in ('witnesses', 'dimensions')}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('propose', 'replay', 'baseline'))
    parser.add_argument('certificate', type=Path)
    parser.add_argument('witness_file', type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    if args.mode == 'propose':
        propose(certificate, args.witness_file)
        return
    started = time.monotonic()
    witnesses = json.loads(args.witness_file.read_text())['witnesses'] if args.mode == 'replay' else None
    result = replay(certificate, psd_witnesses=witnesses)
    seconds = time.monotonic()-started
    previous_path = args.certificate.with_name(args.certificate.stem+'_replay.json')
    previous = json.loads(previous_path.read_text())['lower_replay']
    for field in ('target', 'chain_sites', 'open_lower_density', 'periodic_lower_density', 'local_sum_dimensions'):
        if result[field] != previous[field]:
            raise ValueError('Independent baseline replay differs: '+field)
    groups = ['local_sectors']
    current_blocks = result['local_sectors']+result['overlap']['sectors']+result['joint_overlap']['sectors']
    previous_blocks = previous['local_sectors']+previous['overlap']['sectors']+previous['joint_overlap']['sectors']
    if len(current_blocks) != len(previous_blocks):
        raise ValueError('Block coverage differs')
    for current, old in zip(current_blocks, previous_blocks):
        if current['sector'] != old['sector'] or any(current['psd'][k] != old['psd'][k] for k in ('dimension', 'rank', 'nullity', 'positive_semidefinite')):
            raise ValueError('Exact PSD metadata differs')
    files = {Path(__file__).resolve(), args.certificate.resolve(), previous_path.resolve()}
    if args.mode == 'replay':
        files.add(args.witness_file.resolve())
    for module in tuple(sys.modules.values()):
        source = getattr(module, '__file__', None)
        if source and str(Path(source).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(source).resolve())
    result = {'accepted': True, 'mode': args.mode, 'seconds': seconds, 'lower_replay': result,
              'blocks_checked': len(current_blocks),
              'congruence_blocks': sum('method' in b['psd'] for b in current_blocks),
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
              'scope': 'Fresh standard-library all-block reconstruction and exact PSD acceptance. Matches prior exact lower and all PSD ranks/nullities. Congruence witnesses are optional auxiliary proof data, not a new physical energy bound.'}
    output = args.witness_file.with_name(args.witness_file.stem+'_'+args.mode+'.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('lower_replay', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
