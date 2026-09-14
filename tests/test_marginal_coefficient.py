"""Coefficient-only search contract and exact export checks."""
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
import unittest

from experiments.marginal_coefficient import coefficient_rows,dictionaries,export,solve_coefficients
from experiments.marginal_symbolic import number_shift


class CoefficientTests(unittest.TestCase):
    def test_hermitian_coefficient_dimensions(self):
        self.assertEqual(len(coefficient_rows(6)),352)
        self.assertEqual(len(coefficient_rows(8)),2039)

    def test_direct_number_problem_exports_sound_bound(self):
        h=number_shift(4,0)
        blocks=dictionaries(4,'quadratic')
        solution=solve_coefficients(h,4,2,blocks)
        _,receipt=export(h,4,2,blocks,solution)
        self.assertGreater(F(receipt['lower']),F(1999,1000))
        self.assertLessEqual(F(receipt['lower']),2)

    def test_search_runs_with_old_sector_modules_forbidden(self):
        script="""
import builtins
original=builtins.__import__
def guarded(name,*args,**kwargs):
    if name.startswith('experiments.marginal_hopping'):
        raise AssertionError('Sector dependency forbidden')
    return original(name,*args,**kwargs)
builtins.__import__=guarded
from experiments.marginal_coefficient import dictionaries,solve_coefficients,export
from experiments.marginal_symbolic import number_shift
h=number_shift(4,0); blocks=dictionaries(4,'quadratic')
cert,result=export(h,4,2,blocks,solve_coefficients(h,4,2,blocks))
assert result['lower_float']>1.999
"""
        result=subprocess.run([sys.executable,'-c',script],cwd=Path(__file__).resolve().parents[1],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_invalid_problem_rejected(self):
        with self.assertRaises(ValueError):dictionaries(5,'mixed')
        with self.assertRaises(ValueError):dictionaries(6,'unknown')
        with self.assertRaises(ValueError):solve_coefficients(number_shift(4,0),4,5,dictionaries(4,'quadratic'))


if __name__=='__main__':unittest.main()
