import unittest
from experiments.marginal_symbolic import decode, canonical
from experiments.marginal_transfer_verify import apply_word
from research.global_response_20260913.exceptional.local_coupling import local_blocks
import json

class LocalCouplingTests(unittest.TestCase):
    def test_independent_car_coefficients_all_sources(self):
        for name in ('h6','h8'):
            data=json.load(open(f'results/molecular_collective_20260913/campaign/{name}/fixture.json'))
            got,inputs,outputs=local_blocks(data); m=data['modes']; r=m-4; expected=[[{} for _ in inputs] for _ in outputs]
            for word,co in decode(data['hamiltonian'],m,4).items():
                rest=tuple((c,i) for c,i in word if i<r)
                local=tuple((c,i-r) for c,i in word if i>=r)
                crossings=sum(i>=r and j<r for a,(_,i) in enumerate(word) for _,j in word[a+1:])
                for col,src in enumerate(inputs):
                    v=apply_word(local,src)
                    if v and v[0] in outputs:
                        row=outputs.index(v[0]); parity=crossings+len(local)*(data['particles']-src.bit_count())
                        sign=-1 if parity%2 else 1
                        expected[row][col][rest]=expected[row][col].get(rest,0)+co*v[1]*sign
            self.assertEqual(got,[[canonical(p) for p in row] for row in expected])
            self.assertTrue(any(p for row in got for p in row))

if __name__=='__main__': unittest.main()
