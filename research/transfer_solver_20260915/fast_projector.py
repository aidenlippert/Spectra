"""Measure an existing exact projection algorithm on new molecular row sets."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from research.collective_completion_20260914 import spin_rows
from research.sector_quotient_20260914.fast_twirl import twirl
from research.transfer_solver_20260915.budget import dump


def run(case, output):
    started = time.monotonic()
    metadata = json.loads((case/'prepared/frame.json').read_text())
    rows = [tuple(map(tuple, w)) for w in metadata['rows']]
    spin_rows.twirl = twirl
    T, selected, receipt = spin_rows.build(rows)
    elapsed = time.monotonic()-started
    old = sparse.load_npz(case/'prepared/twirl.npz')
    difference = T-old
    difference.eliminate_zeros()
    same_rows = np.array_equal(selected, np.load(case/'prepared/selected.npy'))
    if difference.nnz or not same_rows:
        raise ValueError('Optimized projector differs from the frozen constructor')
    dump(output, {'seconds': elapsed, 'projector_entries_identical': True,
                   'selected_rows_identical': True, 'exact_integer_projection_checks': receipt,
                   'original_accepting_projector_unchanged': True,
                   'cold_full_pipeline_speedup_measured': False})
    print(json.dumps({'seconds': elapsed, 'identical': True}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    run(a.case.resolve(), a.output.resolve())
