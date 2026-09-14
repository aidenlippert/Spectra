"""Refine an existing connected component certificate with second-order response."""
from fractions import Fraction as F
import json
from pathlib import Path
import time

from experiments.marginal_molecular_gershgorin import matrix
from experiments.marginal_molecular_components import prepare, replay, suggest_lower, strongest_edge_partition
from experiments.marginal_implicit_certificate import rational_text


def refine(path, automatic=False):
    source = Path(path)
    out = source.parent.with_name(source.parent.name + ('_automatic' if automatic else '') + '_second_order')
    if out.exists():
        raise ValueError('Preserve previous second-order export')
    started = time.monotonic()
    c = json.loads(source.read_text())
    old = replay(c)
    states, a = matrix(c)
    if automatic:
        if c.get('component_bound', 'gershgorin') != 'gershgorin':
            raise ValueError('Automatic refinement requires fresh Gershgorin block bounds')
        q = [i for i, state in enumerate(states) if state not in c['retained_states']]
        partition = strongest_edge_partition(a, q, max_size=7)
        c['reference_components'] = [[states[i] for i in group] for group in partition]
        c['reference_selection'] = {'method': 'strongest_edges', 'max_size': 7}
    else:
        partition = [[states.index(s) for s in group] for group in c['reference_components']]
    data = prepare(a, [states.index(s) for s in c['retained_states']], c.get('component_lowers'), partition)
    c['lower'] = rational_text(suggest_lower(data, F(old['upper']), F(old['lower']), 'second_order'))
    c['response'] = 'second_order'
    receipt = replay(c)
    receipt.update(previous_width=old['width'], previous_width_float=old['width_float'],
                   source=str(source), elapsed_seconds=time.monotonic() - started)
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(c, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('width_float', 'previous_width_float', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate')
    parser.add_argument('--automatic', action='store_true')
    args = parser.parse_args()
    refine(args.certificate, args.automatic)
