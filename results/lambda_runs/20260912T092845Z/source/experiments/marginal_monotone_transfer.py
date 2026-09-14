"""Exact transfer of a weighted-row gap under off-diagonal contraction.

The reference must be a checked joint-polynomial metric certificate. This is
not valid for an arbitrary PSD reference proof. The actual target Hamiltonian
and every changed coherent amplitude are reconstructed before acceptance.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from experiments.marginal_polynomial_metric import JointPolynomial, replay as replay_reference


def replay(certificate):
    if certificate.get('kind') != 'joint_polynomial_monotone_transfer_v1':
        raise ValueError('Unsupported monotone-transfer certificate')
    reference = certificate.get('reference_certificate')
    if type(reference) is not dict or reference.get('kind') != 'joint_polynomial_metric_gap_v1':
        raise ValueError('A joint-polynomial weighted-row reference is required')
    if any(certificate.get(key) != reference.get(key) for key in ('modes', 'particles')):
        raise ValueError('Reference and target sectors must agree')
    if type(certificate.get('target_lower')) is not str:
        raise ValueError('Exact target lower bound required')
    gamma = F(certificate['target_lower'])
    if gamma > F(reference['target_lower']):
        raise ValueError('Target exceeds reference gap')
    old = JointPolynomial(reference)
    target = JointPolynomial(dict(certificate, polynomial_metric=reference['polynomial_metric']))
    a, b = old.o, target.o
    if a.oracle.diagonal != b.oracle.diagonal:
        raise ValueError('Monotone transfer requires an unchanged diagonal')
    changes = []
    for key in sorted(a.groups.keys() | b.groups.keys()):
        before, after = a.groups.get(key, {}), b.groups.get(key, {})
        delta = {m: after.get(m, F(0))-before.get(m, F(0))
                 for m in before.keys() | after.keys()}
        delta = {m: v for m, v in delta.items() if v}
        if not delta:
            continue
        if set(delta) != {0}:
            raise ValueError('Only constant grouped-amplitude shifts are supported')
        create, annihilate = key
        closed = a.close(create | annihilate, annihilate)
        if closed is None:
            continue
        lo, hi = a.amplitude_range(before, *closed)
        shift = delta[0]
        endpoint = hi if shift > 0 else lo
        # max_A[(A+shift)^2-A^2] over the checked amplitude enclosure.
        difference = shift*(2*endpoint+shift)
        if difference > 0:
            raise ValueError('Changed amplitude is not uniformly contractive')
        changes.append({'create': create, 'annihilate': annihilate,
                        'reference_range': [str(lo), str(hi)], 'shift': str(shift),
                        'maximum_squared_magnitude_change': str(difference)})
    proof = replay_reference(reference)
    return {'complement_lower': str(gamma), 'reference_gap': proof,
            'changed_groups': changes,
            'determinant_actions': len(a.oracle.cache)+len(b.oracle.cache),
            'target_amplitude_sign_required': False,
            'certifies_reference_target_segment': True,
            'scope': 'Exact off-diagonal contraction transfers this checked weighted-Q-row reference gap with the same positive metric and full valence complement. Diagonal unchanged; constant coherent-amplitude shifts only. Every convex interpolation between reference and target inherits the bound. No new numerator discovery, physical state loop or target fixed-sign assumption.'}


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--verify', type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(replay(json.loads(args.verify.read_text())), indent=2))
