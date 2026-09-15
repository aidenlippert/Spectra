"""Connected hopping stress of a fixed finite molecular block reference."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import time

from experiments.marginal_molecular_gershgorin import matrix, upper, row_bounds
from experiments.marginal_molecular_components import components, prepare, suggest_lower, replay
from experiments.marginal_symbolic import decode, encode, add, adj, mono
from experiments.marginal_symmetry_transfer import integer_upper

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/marginal_molecular_gershgorin/square_coupling/certificate.json'


def perturb(certificate, strength):
    strength = F(strength)
    if not 0 < strength <= F(1, 1000):
        raise ValueError('This stress sweep requires 0 < strength <= 1/1000')
    result = copy.deepcopy(certificate)
    h = decode(result['hamiltonian'], 8, 4)
    for i in range(7):
        term = mono(((1, i), (0, i + 1)), strength)
        h = add(h, term, adj(term))
    result['hamiltonian'] = encode(h)
    return result


def run(strength=F(1, 1000)):
    source = json.loads(SOURCE.read_text())
    c = perturb(source, strength)
    out = ROOT / 'results/marginal_molecular_gershgorin' / ('square_connected_' + str(strength).replace('/', '_'))
    if out.exists():
        raise ValueError('Preserve previous connected stress export')
    started = time.monotonic()
    states, original = matrix(source)
    retained = [states.index(s) for s in source['retained_states']]
    q = [i for i in range(70) if i not in retained]
    partition = components(original, q)
    _, a = matrix(c)
    c['independent_upper'] = integer_upper(decode(c['hamiltonian'], 8, 4), 8, 4)
    u = upper(a, c['independent_upper'])
    data = prepare(a, retained, partition=partition)
    start_lower = min(row_bounds(a, list(range(70))).values()) - 2 * data['off_block_norm'] - 1
    lower = suggest_lower(data, u, start_lower, 'moments')
    c.update(kind='molecular_component_schur_v1', lower=str(lower), component_bound='gershgorin', response='moments',
             reference_components=[[states[i] for i in group] for group in partition])
    receipt = replay(c)
    receipt.update(strength=str(strength), source=str(SOURCE), elapsed_seconds=time.monotonic() - started,
                   actual_q_component_dimensions=[len(g) for g in components(a, q)],
                   upper_method='Full 70-state numerical eigenvector proposed, integer Rayleigh quotient checked exactly')
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(c, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    # Control: retain the actual connected Q inverse, without the block surrogate.
    control = copy.deepcopy(c)
    control.pop('reference_components')
    control['response'] = 'exact'
    exact_data = prepare(a, retained)
    control['lower'] = str(suggest_lower(exact_data, u, lower, True))
    control_receipt = replay(control)
    (out / 'control_certificate.json').write_text(json.dumps(control, indent=2) + '\n')
    (out / 'control_receipt.json').write_text(json.dumps(control_receipt, indent=2) + '\n')
    print(json.dumps({'strength': str(strength), 'width_float': receipt['width_float'],
                      'off_block_norm': receipt['off_block_norm_bound'],
                      'response_degrees': [x['degree'] for x in receipt['response_recurrences']],
                      'control_width': control_receipt['width_float'], 'elapsed_seconds': time.monotonic() - started}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strength', default='1/1000')
    args = parser.parse_args()
    run(F(args.strength))
