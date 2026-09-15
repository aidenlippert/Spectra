"""Check the imported representation against the actual local CAR maps."""
from collections import Counter
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys
import time
import unittest

from research.h8_spin_import_20260914.budget import ROOT, OUT

PACKAGE = OUT/'imported/Spectra_H8_0767448mHa'
if (OUT/'manifest.json').exists():
    raise RuntimeError('Sealed validation pass')
os.environ['SPECTRA_ROOT'] = str(ROOT)
os.environ['SPECTRA_OUTPUT'] = str(OUT/'representation_check')
sys.path.insert(0, str(PACKAGE/'discovery'))
from spin_irrep import raising, kernel, construct, load_operator, load_source
import numpy as np
from experiments.marginal_symbolic import canonical, mono, add, product, scale
from research.collective_completion_20260914.spin_screen import check_sector, check

RECEIPT = {}


class RepresentationChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.op = load_operator()
        cls.groups = cls.op.meta['groups']
        cls.raising, cls.next = raising(cls.groups)

    def test_every_mixed_word_against_local_CAR(self):
        S = add(*(mono(((1, i), (0, i+1))) for i in range(0, 16, 2)))
        checked = 0
        for g in range(42, 58):
            for j, word in enumerate(self.groups[g]['words']):
                B = canonical(mono(tuple(map(tuple, word))))
                literal = add(product(S, B), scale(product(B, S), -1))
                proposed = {}
                if g in self.raising:
                    column = self.raising[g].getcol(j).tocoo()
                    h = self.next[g]
                    for i, value in zip(column.row, column.data):
                        term = canonical(mono(tuple(map(tuple, self.groups[h]['words'][i]))))
                        proposed = add(proposed, scale(term, int(value)))
                self.assertEqual(literal, proposed, (g, j))
                checked += 1
        self.assertEqual(checked, 3840)
        RECEIPT['mixed_words_checked_by_exact_local_CAR'] = checked

    def test_integer_highest_weight_kernel_and_rank(self):
        cases = []
        for g, U in self.raising.items():
            if U.shape != (112, 368):
                continue
            self.assertTrue(np.array_equal((U@U.T).toarray(), 3*np.eye(112, dtype=np.int64)))
            K = kernel(U)
            K6 = np.rint(6*K).astype(np.int64)
            self.assertEqual(K6.shape, (368, 256))
            self.assertTrue(np.array_equal(K, K6/6))
            self.assertTrue(np.all(U@K6 == 0))
            gram = K6.T@K6
            self.assertTrue(np.all(gram == np.diag(np.diag(gram))))
            self.assertTrue(np.all(np.diag(gram) > 0))
            cases.append({'physical_group': g, 'integer_kernel_dimension': 256,
                          'squared_norms_of_six_times_basis': sorted(set(map(int, np.diag(gram))))})
        self.assertEqual(len(cases), 4)
        RECEIPT['exact_kernel_checks'] = cases

    def test_integer_spin_commutator(self):
        for g in range(42, 58):
            n = len(self.groups[g]['words'])
            matrix = np.zeros((n, n), dtype=np.int64)
            if g in self.raising:
                matrix -= (self.raising[g].T@self.raising[g]).toarray()
            for h, target in self.next.items():
                if target == g:
                    matrix += (self.raising[h]@self.raising[h].T).toarray()
            word = self.groups[g]['words'][0]
            two_m = sum((1 if c else -1)*(1 if p % 2 == 0 else -1) for c, p in word)
            self.assertTrue(np.array_equal(matrix, two_m*np.eye(n, dtype=np.int64)))
        RECEIPT['integer_spin_commutator_groups_checked'] = 16

    def test_gram_transport_against_local_coefficient_maps(self):
        op = self.op
        Q, _, _ = load_source(op)
        positions = {members[0]: i for i, members in enumerate(op.members)}
        specs = construct(op, Q)
        reconstructed = sum((op.M[positions[g]]@((V@q@V.T).ravel()) for g, V, q, _ in specs),
                            np.zeros_like(op.rhs))
        expected = op.A(Q)
        error = float(np.max(abs(reconstructed-expected)))
        self.assertLess(error, 1e-12)
        dimensions = [len(q) for _, _, q, _ in specs]
        self.assertEqual(sum(n*n for n in dimensions), 335168)
        RECEIPT.update(warm_start_max_coefficient_error=error, block_dimensions=dimensions,
                       matrix_entries=sum(n*n for n in dimensions),
                       independent_symmetric_entries=sum(n*(n+1)//2 for n in dimensions))
        full = load_operator()
        for i, members in enumerate(full.members):
            if members[0] >= 42:
                n = full.V[i].shape[0]
                full.V[i] = np.eye(n)
                full.identity[i] = True
        rng = np.random.default_rng(20260915)
        errors = []
        for _ in range(2):
            Q = []
            for V in full.V:
                factor = rng.normal(size=(V.shape[1], 2))
                Q.append(factor@factor.T)
            specs = construct(full, Q)
            new = sum((full.M[positions[g]]@((V@q@V.T).ravel()) for g, V, q, _ in specs),
                      np.zeros_like(full.rhs))
            old = full.A(Q)
            error = float(np.linalg.norm(new-old)/np.linalg.norm(old))
            self.assertLess(error, 2e-13)
            errors.append(error)
        RECEIPT['random_full_Gram_transport_relative_errors'] = errors
        RECEIPT['transport_check_scope'] = 'Floating diagnostics on local coefficient maps; energy acceptance is independent exact CAR replay.'

    def test_pruning_preserves_every_other_coefficient(self):
        read = lambda name: json.loads((PACKAGE/'certificates'/name).read_text())
        old, new = read('unpruned_singlet.json'), read('singlet.json')
        old_blocks = old['core'].pop('blocks')
        new_blocks = new['core'].pop('blocks')
        self.assertEqual(old, new)
        self.assertEqual(len(old_blocks), len(new_blocks))
        removed = 0
        bound_numerator = 0
        kept = 0
        for a, b in zip(old_blocks, new_blocks):
            af, bf = a.pop('factor'), b.pop('factor')
            self.assertEqual(a, b)
            ac, bc = Counter(map(tuple, af)), Counter(map(tuple, bf))
            self.assertFalse(bc-ac)
            for row, count in (ac-bc).items():
                norm1 = sum(abs(v) for v in row)
                self.assertLess(norm1, 1000)
                removed += count
                bound_numerator += count*norm1**2
            kept += len(bf)
        bound = F(bound_numerator, old['core']['denominator']**2)
        self.assertEqual((removed, kept), (377, 203))
        self.assertEqual(bound, F(198941, 125000000000000000))
        RECEIPT['pruning'] = {'removed_rows': removed, 'retained_rows': kept,
                              'exact_removed_square_norm_bound_Ha': str(bound),
                              'all_other_certificate_content_identical': True}

    def test_original_checker_rejects_changed_inputs_and_missing_spin_piece(self):
        read = lambda path: json.loads(path.read_text())
        data = read(PACKAGE/'inputs/fixture.json')
        certificate = read(PACKAGE/'certificates/singlet.json')
        wrong = json.loads(json.dumps(certificate))
        wrong['particles'] = 7
        with self.assertRaises(ValueError):
            check_sector(data, wrong)
        wrong = json.loads(json.dumps(certificate))
        wrong['hamiltonian'][0]['coefficient'] = '123'
        with self.assertRaises(ValueError):
            check_sector(data, wrong)
        with self.assertRaises(ValueError):
            check(data, certificate, {})
        RECEIPT['refusals_checked'] = ['changed particle number', 'changed Hamiltonian coefficient', 'missing nonsinglet proof']


if __name__ == '__main__':
    start = time.monotonic()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RepresentationChecks)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    RECEIPT.update(tests_run=result.testsRun, successful=result.wasSuccessful(), seconds=time.monotonic()-start,
                   no_new_numerical_optimization_run=True)
    path = OUT/'representation_check/receipt.json'
    with path.open('x') as stream:
        json.dump(RECEIPT, stream, indent=2)
        stream.write('\n')
    print(json.dumps(RECEIPT, indent=2))
    if not result.wasSuccessful():
        raise SystemExit(1)
