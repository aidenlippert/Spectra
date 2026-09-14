"""Matched finite ablations and a charged full-interval geometry transfer."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys
import time

from research.compact_response_20260913 import program, closure, trial
from research.certificate_scaling.streaming_reference_upper import upper
from research.molecular_collective_20260913.core import digest

OUT = program.ROOT/'results/composable_response_20260913'


def rational_factor(matrix, denominator, eta=F(0)):
    import numpy as np
    d = len(matrix); fden = 10**11
    L = np.linalg.cholesky(np.array(matrix, dtype=float)/float(denominator)-(float(eta)+1e-7)*np.eye(d))
    factor = [[int(round(float(L[i, j])*fden)) for j in range(i+1)] for i in range(d)]
    margin = closure.check_factor(matrix, denominator, eta, factor, fden)
    return factor, str(margin)


def ablations():
    import numpy as np
    start = time.monotonic(); data, tail, reference = program.load_case('h6')
    response = json.loads((program.OUT/'h6_program.json').read_text()); target = F(response['target_Ha']); delta = F(response['delta_Ha'])
    groups, hden, counts = closure.blocks(data); xden = 10**10; records = []
    for label in ('zero_response', 'constant_response', 'polynomial_response', 'direct_full_factor'):
        t = time.monotonic(); blocks = []; failures = []; actions = 0
        for g in groups:
            A, Bt, D, den = closure.shifted(g, hden, target); n = len(A); q = len(D)
            if label == 'direct_full_factor':
                full_den = den; multiplier = den//hden
                matrix = [[h*multiplier-(int(target*den) if i == j else 0) for j, h in enumerate(row)] for i, row in enumerate(g['H'])]
                factor, margin = rational_factor(matrix, full_den)
                blocks.append({'key': g['key'], 'factor': factor, 'margin_Ha': margin}); continue
            bf = np.array(Bt, dtype=float).reshape(q, n)/float(den); df = np.array(D, dtype=float).reshape(q, q)/float(den)
            if label == 'zero_response':
                xf = np.zeros((q, n))
            elif label == 'constant_response':
                xf = 2*bf/float(F(response['M_Ha'])+delta)
            else:
                xf = program.apply_response(lambda v: df@v, bf, response) if q else np.zeros((0, n))
                actions += (response['order']-1)*n if q else 0
            X = np.rint(xf*xden).astype(np.int64).tolist()
            if label == 'constant_response':
                X = [[round(F(2*c*xden, den)/(F(response['M_Ha'])+delta)) for c in row] for row in Bt]
            K, kden, actual = closure.residual_and_K(A, Bt, D, den, X, xden, delta)
            eta = program.ceil_grid(actual, 10**12)
            M = np.array(K, dtype=float)/float(kden)-float(eta)*np.eye(n)
            e, U = np.linalg.eigh(M)
            if e[0] <= 1e-7:
                v = np.rint(U[:, 0]*10**9).astype(np.int64).tolist()
                value = F(sum(v[i]*K[i][j]*v[j] for i in range(n) for j in range(n)), kden)-eta*sum(x*x for x in v)
                assert value < 0
                failures.append({'key': g['key'], 'vector': v, 'X': X, 'eta_Ha': str(eta), 'exact_negative_expectation_Ha': str(value/sum(x*x for x in v))})
            else:
                factor, margin = rational_factor(K, kden, eta)
                blocks.append({'key': g['key'], 'X': X, 'eta_Ha': str(eta), 'factor': factor, 'margin_Ha': margin})
        record = {'kind': label, 'fixture_sha256': digest(data), 'response_program_sha256': digest(response),
            'target_Ha': str(target), 'response_denominator': xden, 'factor_denominator': 10**11,
            'blocks': blocks, 'failures': failures, 'all_blocks_pass': not failures,
            'D_actions_on_columns': actions, 'seconds': time.monotonic()-t}
        raw = json.dumps(record, separators=(',', ':'))+'\n'; path = OUT/f'{label}.json'; path.write_text(raw)
        records.append({'control': label, 'accepted': not failures, 'failed_blocks': len(failures),
            'worst_negative_expectation_Ha': min((float(F(x['exact_negative_expectation_Ha'])) for x in failures), default=None),
            'certificate_or_obstruction_bytes': path.stat().st_size, 'seconds': record['seconds'], 'D_actions_on_columns': actions})
    return {'controls': records, 'shared_block_construction': counts, 'wall_seconds': time.monotonic()-start}


def transfer():
    import numpy as np
    start = time.monotonic(); data, tail, hf = trial.load_case('fresh_h6_1p6')
    groups, hden, counts = closure.blocks(data)
    group = next(g for g in groups if 63 in g['states'])
    values, vectors = np.linalg.eigh(np.array(group['H'], dtype=float)/float(hden))
    amps = np.rint(vectors[:, 0]*10**10).astype(np.int64).tolist()
    witness = {'states': [s for s, a in zip(group['states'], amps) if a], 'amplitudes': [a for a in amps if a]}
    U, upper_receipt = upper(data, witness); reference = {'independent_upper': witness, 'upper': str(U)}
    (OUT/'transfer_upper.json').write_text(json.dumps(reference, indent=2)+'\n')
    response, rstats = program.propose(data, tail, U-F(14, 10000))
    (OUT/'transfer_program.json').write_text(json.dumps(response, indent=2)+'\n')
    cert, obstruction, stats = closure.propose(data, tail, response)
    (OUT/'transfer_closure.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n')
    return {'case': 'fresh_h6_1p6', 'upper_discovery': 'Full diagonalization of the 200-state block containing HF; enumerated FCI reference work.',
        'reference_block_dimension': len(group['states']), 'reference_basis_construction': counts,
        'upper_replay': upper_receipt, 'program_discovery': rstats, 'closure_discovery': stats,
        'new_full_width_mHa': 1.4, 'target_Ha': response['target_Ha'], 'upper_Ha': str(U),
        'certificate_bytes': (OUT/'transfer_closure.json').stat().st_size,
        'wall_seconds': time.monotonic()-start, 'scope': 'Full-interval transfer using the unchanged response rule and explicitly enumerated upper and retained proof.'}


def replay():
    start = time.monotonic(); data, tail, reference = program.load_case('h6')
    response = json.loads((program.OUT/'h6_program.json').read_text()); program.check(data, tail, response)
    groups, hden, counts = closure.blocks(data); rows = []
    for label in ('zero_response', 'constant_response', 'polynomial_response', 'direct_full_factor'):
        cert = json.loads((OUT/f'{label}.json').read_text())
        if cert.get('kind') != label or any(type(cert[k]) is not int or not 0 < cert[k] <= 10**14 for k in ('response_denominator','factor_denominator')):
            raise ValueError('Invalid matched-control kind or denominator')
        if cert['fixture_sha256'] != digest(data) or cert['response_program_sha256'] != digest(response) or F(cert['target_Ha']) != F(response['target_Ha']):
            raise ValueError('Matched-control input binding failed')
        if any('vector' in row or 'factor' not in row for row in cert['blocks']) or any('vector' not in row or 'factor' in row for row in cert['failures']):
            raise ValueError('Factor successes and negative witnesses must be classified correctly')
        if label == 'direct_full_factor' and cert['failures']:
            raise ValueError('Direct control requires complete factors')
        supplied = cert['blocks']+cert['failures']
        if sorted(r['key'] for r in supplied) != sorted(g['key'] for g in groups):
            raise ValueError('Matched control did not cover all conserved blocks')
        for row in supplied:
            g = next(g for g in groups if g['key'] == row['key'])
            A, Bt, D, den = closure.shifted(g, hden, F(response['target_Ha']))
            if label == 'direct_full_factor':
                K = [[h*(den//hden)-(int(F(response['target_Ha'])*den) if i == j else 0) for j, h in enumerate(r)] for i, r in enumerate(g['H'])]
                kden = den; eta = F(0)
            else:
                K, kden, actual_eta = closure.residual_and_K(A, Bt, D, den, row['X'], cert['response_denominator'], F(response['delta_Ha']))
                eta = F(row['eta_Ha'])
                if eta < actual_eta:
                    raise ValueError('Understated residual penalty')
                if label == 'zero_response' and any(x for r in row['X'] for x in r):
                    raise ValueError('Zero-response control contains a response')
                if label == 'constant_response':
                    expected = [[round(F(2*c*cert['response_denominator'], den)/(F(response['M_Ha'])+F(response['delta_Ha']))) for c in r] for r in Bt]
                    if row['X'] != expected:
                        raise ValueError('Constant-response control is not the specified constant')
            if 'vector' in row:
                v = row['vector']
                if len(v) != len(K) or any(type(x) is not int for x in v) or not any(v):
                    raise ValueError('Invalid control obstruction')
                value = F(sum(v[i]*K[i][j]*v[j] for i in range(len(v)) for j in range(len(v))), kden)/sum(x*x for x in v)-eta
                if value >= 0 or str(value) != row['exact_negative_expectation_Ha']:
                    raise ValueError('Control obstruction failed')
            elif str(closure.check_factor(K, kden, eta, row['factor'], cert['factor_denominator'])) != row['margin_Ha']:
                raise ValueError('Control factor margin changed')
        if cert['all_blocks_pass'] != (not cert['failures']):
            raise ValueError('Invalid matched-control outcome')
        rows.append({'control': label, 'accepted': cert['all_blocks_pass'], 'failed_blocks': len(cert['failures'])})
    data, tail, _ = trial.load_case('fresh_h6_1p6')
    transfer = closure.check(data, tail, json.loads((OUT/'transfer_program.json').read_text()),
        json.loads((OUT/'transfer_closure.json').read_text()), json.loads((OUT/'transfer_upper.json').read_text()))
    return {'controls': rows, 'ablation_shared_block_construction': counts, 'transfer': transfer,
        'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': [x for x in ('numpy', 'scipy', 'cvxpy', 'pyscf') if x in sys.modules]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true'); args = parser.parse_args()
    if args.replay:
        result = replay(); assert not result['numerical_packages_loaded']; name = 'finite_replay.json'
    else:
        result = {'ablations': ablations(), 'transfer': transfer()}; name = 'finite_discovery.json'
    (OUT/name).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': name, 'transfer_width_mHa': result['transfer'].get('width_mHa', result['transfer'].get('new_full_width_mHa'))}), flush=True)
