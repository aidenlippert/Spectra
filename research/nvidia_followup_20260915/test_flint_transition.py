"""A warm standard-library CAR cache must not leak into a FLINT lower replay."""
import subprocess
import sys
import unittest


class TransitionTests(unittest.TestCase):
    def test_warm_upper_style_cache_then_lower_arithmetic(self):
        code='''
from fractions import Fraction
from experiments import marginal_symbolic as symbolic
from research.collective_completion_20260914.spin_screen import spin_squared
before={w:str(v) for w,v in spin_squared(6).items()}
symbolic.product({((1,0),(0,0)):Fraction(2,3)},{():Fraction(7,11)})
assert symbolic.word_product.cache_info().currsize > 0
from research.nvidia_followup_20260915.replay import enable_flint_rationals
enable_flint_rationals()
assert symbolic.word_product.cache_info().currsize == 0
after={w:str(v) for w,v in spin_squared(6).items()}
assert before==after
import flint
got=symbolic.product({((1,0),(0,0)):flint.fmpq(2,3)},{():flint.fmpq(7,11)})
assert str(next(iter(got.values())))=='14/33'
assert 'numpy' not in __import__('sys').modules
'''
        result=subprocess.run([sys.executable,'-B','-c',code],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)


if __name__=='__main__':
    unittest.main()
