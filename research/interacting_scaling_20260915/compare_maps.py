"""Require identical coefficient systems when measuring a construction change."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from research.interacting_scaling_20260915.budget import dump


def compare(reference, candidate, output):
    a = json.loads((reference/'frame.json').read_text())
    b = json.loads((candidate/'frame.json').read_text())
    for key in ('fixture_sha256', 'rows', 'groups', 'blocks', 'number_basis', 'spin_basis', 'ladder_basis'):
        if a[key] != b[key]: raise ValueError(('Declared coefficient system changed', key))
    matrices = ['twirl', 'free', 'normal']+[f'map_{i}' for i in range(len(a['blocks']))]
    for name in matrices:
        left, right = [sparse.load_npz(p/(name+'.npz')) for p in (reference, candidate)]
        difference = left-right; difference.eliminate_zeros()
        if difference.nnz: raise ValueError(('Coefficient matrix changed', name, float(np.max(abs(difference.data)))))
    for name in ('rhs', 'scale', 'weights', 'selected', 'moments', 'physical_dual'):
        if not np.array_equal(np.load(reference/(name+'.npy')), np.load(candidate/(name+'.npy'))):
            raise ValueError(('Numerical input changed', name))
    dump(output, {'all_maps_and_inputs_identical': True, 'sparse_matrices_compared': len(matrices),
        'reference_constructor_seconds': a['total_seconds'], 'candidate_constructor_seconds': b['total_seconds'],
        'construction_ratio': a['total_seconds']/b['total_seconds'],
        'scope': 'Matched operator preparation, excluding shared state discovery, optimization and verification'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('reference', type=Path); p.add_argument('candidate', type=Path); p.add_argument('output', type=Path)
    a = p.parse_args(); compare(a.reference, a.candidate, a.output)
