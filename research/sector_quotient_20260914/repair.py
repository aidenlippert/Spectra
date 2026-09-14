"""Exact quotient projection and an optional proposed residual wedge witness.

No remainder is discarded. An independent replay must accept the resulting
certificate; this merely chooses a quartic multiplier from the actual squares.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import verified_residual, decode, encode, add, product, number_shift, scale
from research.sector_quotient_20260914.fast_twirl import twirl
from research.sector_quotient_20260914.algebra import from_polynomial, to_polynomial, lift, combine, split_three_body, inner
from research.reconstruction_compression_20260914.inputs import sha, dump


def run(source, out):
    from research.sector_quotient_20260914.budget import OUT
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed campaign')
    start = time.monotonic(); source = Path(source); folder = Path(out); folder.mkdir(parents=True, exist_ok=False)
    cert = json.loads(source.read_text()); core = cert['core']; m = cert['modes']
    if not cert.get('spin_twirl') or not cert.get('singlet'):
        raise ValueError('Expected the declared spin-averaged singlet certificate')
    residual, rec = verified_residual(core)
    print(json.dumps({'stage': 'expanded_actual_factors', 'seconds': time.monotonic()-start}), flush=True)
    averaged = twirl(residual)
    X = decode(core['number_multiplier'], m, 4)
    X4 = from_polynomial(twirl({w: c for w, c in X.items() if len(w) == 4}), m, 2)
    R6 = from_polynomial(averaged, m, 3)
    W = combine((-1, R6), (-1, lift(X4, m, 2)))
    V, Z = split_three_body(W, m)
    correction = combine((-1, V), (-1, X4))
    core['number_multiplier'] = encode(add(X, to_polynomial(correction, m, 2)))
    # Exact linear update of the already validated residual. The independent
    # accepting process still re-expands the modified certificate from scratch.
    modified = add(averaged, scale(product(number_shift(m, cert['particles']), to_polynomial(correction, m, 2)), -1))
    if from_polynomial(modified, m, 3) != combine((-1, Z)):
        raise AssertionError('Projected multiplier did not retain exactly -Z')
    from research.certificate_scaling.wedge_spectral_bound import propose
    cert['residual_wedge_witness'] = propose(modified, m, cert['particles'])
    with (folder/'certificate.json').open('x') as f:
        json.dump(cert, f, separators=(',', ':'))
    record = {'source': str(source), 'source_sha256': sha(source), 'exact_quotient_projection': True,
              'exact_sextic_residual_equals_negative_Z': True, 'Z_squared_Frobenius_norm': str(inner(Z, Z)),
              'W_squared_Frobenius_norm': str(inner(W, W)),
              'V_squared_Frobenius_norm': str(inner(V, V)),
              'W_nonzero_coefficients': len(W), 'V_nonzero_coefficients': len(V),
              'Z_nonzero_coefficients': len(Z),
              'Z_is_zero': not Z, 'old_averaged_residual_l1': str(sum(abs(c) for c in averaged.values())),
              'new_averaged_residual_l1': str(sum(abs(c) for c in modified.values())),
              'number_multiplier_correction_norm_squared': str(inner(correction, correction)),
              'many_body_states_enumerated': 0, 'seconds': time.monotonic()-start,
              'proposal_projection': 'fast exact canonical spin derivation, tested against the sealed projector',
              'status': 'requires_independent_exact_replay'}
    dump(folder/'repair.json', record); print(json.dumps(record), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('source'); p.add_argument('out'); a = p.parse_args(); run(a.source, a.out)
