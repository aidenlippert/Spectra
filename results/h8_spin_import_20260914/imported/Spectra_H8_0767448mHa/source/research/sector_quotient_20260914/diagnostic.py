"""Exact diagnosis of an inherited full cubic proof, explicitly not training."""
import json
from fractions import Fraction as F
from math import sqrt
from pathlib import Path
import time
from experiments.marginal_symbolic import verified_residual, decode, product, number_shift, add, scale, mono, encode
from research.sector_quotient_20260914.algebra import from_polynomial, to_polynomial, lift, split_three_body, combine, inner
from research.sector_quotient_20260914.budget import ROOT, OUT
from research.reconstruction_compression_20260914.inputs import sha, dump


def run():
    start = time.monotonic()
    frozen = json.loads((OUT/'frozen_inputs.json').read_text())['h8']
    source = ROOT/frozen['teacher']
    if sha(source) != frozen['sha256']['teacher']:
        raise ValueError('Strong certificate changed')
    core = json.loads(source.read_text())
    data = json.loads((ROOT/frozen['fixture']).read_text())
    m, n = core['modes'], core['particles']
    h = decode(core['hamiltonian'], m, 4)
    if h != decode(data['hamiltonian'], m, 4):
        raise ValueError('Strong certificate Hamiltonian mismatch')
    residual, receipt = verified_residual(core)
    X = decode(core['number_multiplier'], m, 4)
    ideal = product(number_shift(m, n), X)
    S = add(h, mono((), -F(core['b'])), scale(ideal, -1), scale(residual, -1))
    W = from_polynomial(S, m, 3)
    X4 = from_polynomial(X, m, 2)
    R6 = from_polynomial(residual, m, 3)
    mismatch = combine((1, W), (1, lift(X4, m, 2)))
    if mismatch != combine((-1, R6)):
        raise AssertionError('Actual sextic residual does not match W+L(X4)')
    V, Z = split_three_body(W, m)
    VR, ZR = split_three_body(R6, m)
    if combine((1, V), (1, X4), (1, VR)) or combine((1, Z), (1, ZR)):
        raise AssertionError('Projected strong-proof identity failed')
    from research.certificate_scaling.wedge_spectral_bound import replay_residual
    witness = json.loads((ROOT/frozen['teacher_residual']).read_text())
    wedge = replay_residual(core, residual, witness)
    from research.collective_completion_20260914.spin_replay import setup
    _, _, delta = setup(data)
    extracted = {'W_sextic': encode(to_polynomial(W, m, 3)), 'projected_V_quartic': encode(to_polynomial(V, m, 2)),
                 'contraction_free_Z_sextic': encode(to_polynomial(Z, m, 3)), 'old_X4': encode(to_polynomial(X4, m, 2))}
    dump(OUT/'strong_h8_extraction.json', extracted)
    stats = lambda A: {'nonzeros': len(A), 'squared_Frobenius_norm': str(inner(A, A)),
                       'Frobenius_norm_float': sqrt(float(inner(A, A))),
                       'coefficient_l1': str(sum(abs(v) for v in A.values()))}
    rec = {'source': str(source), 'source_sha256': sha(source), 'role': 'diagnostic only',
           'modes': m, 'particles': n, 'quartic_multiplier_terms': len(X4),
           'W': stats(W), 'V': stats(V), 'Z': stats(Z), 'V_plus_X4': stats(combine((1, V), (1, X4))),
           'W_plus_lift_X4': stats(mismatch),
           'exact_W_plus_lift_X4_equals_negative_R6': True,
           'exact_Z_equals_negative_projected_R6': True,
           'exact_lift_membership': not Z,
           'interpretation': 'The nonzero numerical-export remainder is retained and paid by the existing residual certificate.',
           'original_replay': receipt, 'strong_wedge_replay': wedge,
           'strong_certificate_spin_scope': 'Original H on the entire fixed-N sector; no spin averaging or separate screen required.',
           'new_screened_constructor_original_H_spin_defect': str(delta),
           'many_body_states_enumerated': 0, 'seconds': time.monotonic()-start}
    dump(OUT/'strong_h8_diagnostic.json', rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('original_replay', 'strong_wedge_replay')}, indent=2))


if __name__ == '__main__':
    run()
