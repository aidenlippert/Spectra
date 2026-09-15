"""Hamiltonian-only selected CI with sparse Ritz solves and exact upper replay.

Candidate selection is heuristic. Every delivered upper endpoint is an exact
Rayleigh quotient of an explicitly stored integer vector for the original H.
No lower endpoint is discovered here and no asymptotic accuracy is claimed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import lcm
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from experiments.marginal_determinant_tree import DeterminantOracle
from research.certificate_scaling.streaming_reference_upper import compile_term, upper


class CompiledHamiltonian:
    def __init__(self, fixture):
        original = DeterminantOracle(fixture)
        self.original = original
        self.denominator = lcm(*(v.denominator for v in original.h.values()))
        self.terms = [compile_term(w, int(v*self.denominator)) for w, v in original.h.items()]
        self.terms = [t for t in self.terms if t is not None]
        self.diagonal_terms = [t for t in self.terms if t[2] == 0]
        self.cache = {}
        self.term_state_checks = 0
        self.diagonal_checks = 0

    def action(self, state):
        if state not in self.cache:
            row = {}
            for required, occupied, flip, parity, coefficient in self.terms:
                if state & required == occupied:
                    target = state ^ flip
                    signed = -coefficient if (state & parity).bit_count() % 2 else coefficient
                    row[target] = row.get(target, 0) + signed
            self.term_state_checks += len(self.terms)
            self.cache[state] = {s: c/self.denominator for s, c in row.items() if c}
        return self.cache[state]

    def diagonal(self, state):
        result = 0
        for required, occupied, _, parity, coefficient in self.diagonal_terms:
            if state & required == occupied:
                result += -coefficient if (state & parity).bit_count() % 2 else coefficient
        self.diagonal_checks += len(self.diagonal_terms)
        return result / self.denominator


def restricted_matrix(compiled, basis):
    import numpy as np
    from scipy.sparse import coo_matrix
    index = {s: i for i, s in enumerate(basis)}
    rows, cols, data = [], [], []
    for j, state in enumerate(basis):
        for target, value in compiled.action(state).items():
            i = index.get(target)
            if i is not None:
                rows.append(i); cols.append(j); data.append(value)
    matrix = coo_matrix((data, (rows, cols)), shape=(len(basis), len(basis))).tocsr()
    asymmetry = matrix - matrix.T
    if asymmetry.nnz and np.max(np.abs(asymmetry.data)) > 1e-12:
        raise ValueError('Projected Hamiltonian is not symmetric')
    return matrix


def ritz(matrix, previous=None):
    import numpy as np
    from scipy.sparse.linalg import eigsh
    if matrix.shape[0] <= 4:
        values, vectors = np.linalg.eigh(matrix.toarray())
        return float(values[0]), vectors[:, 0]
    if previous is None:
        previous = np.ones(matrix.shape[0])
    values, vectors = eigsh(matrix, k=1, which='SA', tol=1e-10,
                           maxiter=4000, v0=previous)
    return float(values[0]), vectors[:, 0]


def run(fixture_path, out, budgets, reference_path=None, seconds=480):
    import numpy as np
    if not budgets or any(type(x) is not int or x <= 0 or x > 32768 for x in budgets):
        raise ValueError('Support budgets must be positive integers <=32768')
    if sorted(set(budgets)) != budgets or seconds <= 0:
        raise ValueError('Increasing distinct budgets and positive deadline required')
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    raw = Path(fixture_path).read_bytes(); fixture = json.loads(raw)
    compiled = CompiledHamiltonian(fixture)
    if reference_path:
        reference_bytes = Path(reference_path).read_bytes()
        reference = json.loads(reference_bytes)['independent_upper']
        if len(reference['states']) != 1:
            raise ValueError('Reference must be a single supplied HF determinant')
        # Validate reference with the same physical Hamiltonian before use.
        upper(fixture, reference)
        initial = reference['states'][0]
    else:
        initial = (1 << fixture['particles']) - 1
        reference_bytes = None
    if not compiled.original.valid_state(initial):
        raise ValueError('Reference violates target sector')
    basis = [initial]
    previous_coeff = {initial: 1.0}
    rows = []
    candidate_total = 0
    stop = 'budgets_completed'
    for step, budget in enumerate(budgets):
        if time.monotonic()-started > seconds:
            stop = 'deadline_before_iteration'; break
        if len(basis) < budget:
            residual = {}
            energy = sum(previous_coeff[s]*compiled.diagonal(s)*previous_coeff[s]
                         for s in basis) if not rows else rows[-1]['numeric_ritz']
            present = set(basis)
            # Sum signed contributions before ranking, rather than their magnitudes.
            for state, coefficient in previous_coeff.items():
                if abs(coefficient) < 1e-14:
                    continue
                for target, value in compiled.action(state).items():
                    if target not in present:
                        residual[target] = residual.get(target, 0.0) + value*coefficient
            candidate_total += len(residual)
            scored = [(abs(value)/max(abs(compiled.diagonal(s)-energy), 1e-3), s)
                      for s, value in residual.items() if abs(value) > 1e-14]
            scored.sort(key=lambda x: (-x[0], x[1]))
            basis = sorted(present | {s for _, s in scored[:budget-len(basis)]})
        built_at = time.monotonic()
        matrix = restricted_matrix(compiled, basis)
        assembly_seconds = time.monotonic()-built_at
        initial_vector = np.array([previous_coeff.get(s, 0.0) for s in basis])
        solved_at = time.monotonic()
        energy, vector = ritz(matrix, initial_vector)
        solve_seconds = time.monotonic()-solved_at
        previous_coeff = dict(zip(basis, vector))
        residual_norm = float(np.linalg.norm(matrix@vector-energy*vector))
        amplitudes = [int(round(float(x)*10**12)) for x in vector]
        witness = {'states': [s for s, a in zip(basis, amplitudes) if a],
                   'amplitudes': [a for a in amplitudes if a]}
        exact, check = upper(fixture, witness)
        witness_path = out/f'step_{step:02d}_upper.json'
        witness_raw = (json.dumps({'independent_upper': witness, 'upper': str(exact)},
                                 separators=(',', ':'))+'\n').encode()
        witness_path.write_bytes(witness_raw)
        row = {'step': step, 'requested_support': budget, 'basis_size': len(basis),
               'witness_support': len(witness['states']), 'numeric_ritz': energy,
               'projected_residual_norm': residual_norm, 'upper': str(exact),
               'upper_float': float(exact), 'matrix_nnz': matrix.nnz,
               'dense_equivalent_entries': len(basis)**2,
               'assembly_seconds': assembly_seconds, 'solve_seconds': solve_seconds,
               'upper_replay': check, 'elapsed_seconds': time.monotonic()-started,
               'term_state_checks_cumulative': compiled.term_state_checks,
               'diagonal_term_checks_cumulative': compiled.diagonal_checks,
               'candidate_count_cumulative': candidate_total,
               'cached_action_states': len(compiled.cache),
               'witness_file': witness_path.name, 'witness_sha256': sha256(witness_raw).hexdigest()}
        rows.append(row)
        (out/'progress.json').write_text(json.dumps(rows, indent=2)+'\n')
        print(json.dumps(row), flush=True)
        # A saturated connected component is a finite exact-oracle limit, not scaling.
        if len(basis) < budget and step > 0 and len(basis) == rows[-2]['basis_size']:
            stop = 'connected_support_saturated'; break
    result = {'kind': 'signed_residual_sparse_selected_ci_v1',
              'fixture_path': str(fixture_path), 'fixture_sha256': sha256(raw).hexdigest(),
              'modes': fixture['modes'], 'particles': fixture['particles'],
              'reference_state': initial,
              'reference_sha256': sha256(reference_bytes).hexdigest() if reference_bytes else None,
              'reference_source': 'supplied_single_determinant' if reference_bytes else 'lowest_modes_occupation',
              'hamiltonian_terms': len(compiled.terms), 'budgets': budgets,
              'full_fock_enumeration': False, 'saved_FCI_used': False,
              'proposal_drop_threshold': 1e-14, 'denominator_floor': 1e-3,
              'seconds_budget': seconds, 'stop_reason': stop,
              'wall_seconds': time.monotonic()-started, 'rows': rows,
              'scope': 'Heuristic selected-CI discovery; exact upper replay for original finite H. No scalable accuracy theorem.'}
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--budgets', type=int, nargs='+', default=[64,128,256,512,1024,2048,4096])
    parser.add_argument('--seconds', type=float, default=480)
    args = parser.parse_args()
    run(args.fixture, args.out, args.budgets, args.reference, args.seconds)
