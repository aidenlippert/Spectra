"""Apply the generic exact quotient to frozen operator dictionaries, not energies."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

from experiments.marginal_symbolic import canonical
from research.general_mechanism_20260915.sector_quotient import operator_quotient


def encode(value):
    if isinstance(value, F):
        return str(value)
    raise TypeError(type(value).__name__)


def run(prepared, blocks, output):
    started = time.monotonic()
    if output.exists() or output.with_suffix('.run.json').exists():
        raise ValueError('Refusing to overwrite an earlier receipt')
    output.parent.mkdir(parents=True, exist_ok=True)
    ledger = {'scope': 'Frozen operator dictionary quotient; not a cold molecular calculation',
              'status': 'running', 'prepared': str(prepared), 'blocks_requested': blocks}
    try:
        import numpy as np
        sources = {name: hashlib.sha256((prepared/name).read_bytes()).hexdigest()
                   for name in ('frame.json', 'bases.npz')}
        frame = json.loads((prepared/'frame.json').read_text())
        results = []
        with np.load(prepared/'bases.npz') as bases:
            for index in blocks:
                block = frame['blocks'][index]
                words = frame['groups'][block['physical_group']]['words']
                values = bases[f'V_{index}']
                denominator = block.get('integer_denominator', 1)
                if type(denominator) is not int or denominator < 1 or not np.isfinite(values).all():
                    raise ValueError('Invalid rational basis grid')
                numerators = np.rint(denominator*values)
                if not np.array_equal(numerators/denominator, values):
                    raise ValueError('Basis does not round-trip on its declared rational grid')
                operators = []
                for column in range(values.shape[1]):
                    terms = {}
                    for row in np.flatnonzero(numerators[:, column]):
                        word = tuple(tuple(letter) for letter in words[int(row)])
                        terms[word] = terms.get(word, F(0))+F(int(numerators[row, column]), denominator)
                    operators.append(canonical(terms))
                result = operator_quotient(operators, frame['modes'], frame['particles'])
                result['source_block'] = index
                result['basis_denominator'] = denominator
                result['basis_grid_binary_round_trip_checked'] = True
                results.append(result)
                print(json.dumps({k: v for k, v in result.items() if k != 'groups'}), flush=True)
        receipt = {'scope': 'Exact physical singlet-action quotient of supplied rational dictionaries',
                   'input_sha256': sources, 'results': results,
                   'no_global_coefficient_maps_loaded': True,
                   'no_Hamiltonian_or_MPS_loaded': True,
                   'SOS_ideal_equivalence_or_energy_gain_claimed': False}
        with output.open('x') as stream:
            json.dump(receipt, stream, default=encode, indent=2)
            stream.write('\n')
        ledger['status'] = 'completed'
        ledger['receipt_bytes'] = output.stat().st_size
    except Exception as error:
        ledger['status'] = 'failed'
        ledger['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        ledger['measured_function_seconds'] = time.monotonic()-started
        ledger['peak_process_RSS_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
        ledger['timing_excludes'] = ['Interpreter and pre-run imports', 'Inherited dictionary discovery and preparation', 'Test and development work']
        with output.with_suffix('.run.json').open('x') as stream:
            json.dump(ledger, stream, indent=2)
            stream.write('\n')
        print(json.dumps(ledger), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('prepared', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--blocks', type=int, nargs='+', default=[21, 22, 23, 24])
    args = parser.parse_args()
    run(args.prepared.resolve(), args.blocks, args.output.resolve())
