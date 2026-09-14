"""Exact global operator-norm bounds from two CAR sum-of-squares identities."""
from fractions import Fraction as F
from experiments.marginal_symbolic import add, mono, scale, expand_squares


def verify_norm(delta, modes, certificate):
    if type(certificate) is not dict:
        raise ValueError('Explicit two-sided norm certificate required')
    bound = F(certificate['bound'])
    if bound < 0:
        raise ValueError('Nonnegative operator-norm bound required')
    rows = 0
    for key, sign in (('plus', 1), ('minus', -1)):
        recipe = certificate.get(key)
        if type(recipe) is not list or not 1 <= len(recipe) <= 16:
            raise ValueError('Bounded weighted square recipe required')
        total = {}
        for term in recipe:
            if type(term) is not dict or F(term['weight']) < 0:
                raise ValueError('Nonnegative rational square weights required')
            block = term.get('block')
            if (type(block) is not dict or type(block.get('words')) is not list
                    or not 1 <= len(block['words']) <= 64 or type(block.get('factor')) is not list
                    or not 1 <= len(block['factor']) <= 64):
                raise ValueError('Bounded square factor required')
            value, receipt = expand_squares([block], 1, modes, max_degree=2)
            total = add(total, scale(value, F(term['weight'])))
            rows += receipt['factor_rows']
        if add(mono((), bound), scale(delta, sign), scale(total, -1)):
            raise ValueError('Operator-norm square identity does not match the perturbation')
    return bound, {'square_rows': rows, 'scope': 'Exact identities bound I +/- delta = sum positive squares on the entire Fock space'}


def hopping_norm(strength):
    t = F(strength)
    orientation = 1 if t >= 0 else -1
    c = {'bound': str(2*abs(t))}
    for name, sign in (('plus', orientation), ('minus', -orientation)):
        c[name] = [{'weight': str(abs(t)/2), 'block': {
            'words': [[[creation, i]] for i in (0,2,1,3)],
            'factor': [[1, sign if creation == 0 else -sign, 0, 0],
                       [0, 0, 1, sign if creation == 0 else -sign]]}}
            for creation in (0,1)]
    return c


def replay(certificate):
    from experiments.marginal_determinant_tree import DeterminantOracle
    if certificate.get('kind') != 'saturated_operator_norm_v1':
        raise ValueError('Unsupported saturated norm certificate')
    oracle = DeterminantOracle({'modes': certificate.get('modes'), 'particles': certificate.get('particles'),
                                'hamiltonian': certificate['perturbation']})
    bound, receipt = verify_norm(oracle.h, oracle.modes, certificate.get('norm_certificate'))
    rayleigh = oracle.upper(certificate.get('saturating_witness'))
    if abs(rayleigh) != bound:
        raise ValueError('Witness does not saturate the norm bound')
    return dict(receipt, operator_norm=str(bound), saturating_rayleigh=str(rayleigh),
                witness_action_states=len(oracle.cache),
                scope='Exact global Fock-space norm upper bound and a fixed-N physical witness attaining its absolute value; therefore the fixed-N and global operator norms equal the bound.')


def replay_q(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    if certificate.get('kind') != 'saturated_q_operator_norm_v1':
        raise ValueError('Unsupported compressed operator-norm certificate')
    oracle = SpinZeroOracle({'modes':certificate.get('modes'), 'particles':certificate.get('particles'),
                            'hamiltonian':certificate['perturbation']})
    retained = oracle.retained(certificate.get('retained_states'))
    witness = certificate.get('saturating_witness')
    value = oracle.upper(witness)
    if any(s in retained for s, x in zip(witness['states'],witness['amplitudes']) if x):
        raise ValueError('Compressed norm saturation witness must lie entirely in Q')
    result = replay(dict(certificate,kind='saturated_operator_norm_v1'))
    if abs(value) != F(result['operator_norm']):
        raise ValueError('Spin-sector witness does not saturate the compressed norm')
    result.update(retained_dimension=len(retained), witness_support=sum(bool(x) for x in witness['amplitudes']),
        scope='Exact global norm upper bound and a physical Sz=0 witness entirely in the specified Q attaining it. The Q-compressed, Sz=0, fixed-N, and full-Fock operator norms are therefore equal. No complementary configuration enumeration.')
    return result


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify')
    parser.add_argument('--verify-q')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.verify_q:
        print(json.dumps(replay_q(json.loads(Path(args.verify_q).read_text())), indent=2))
    else:
        parser.error('Specify --verify or --verify-q')
