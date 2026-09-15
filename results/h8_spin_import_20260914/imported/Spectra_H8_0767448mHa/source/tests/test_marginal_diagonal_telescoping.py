import unittest
from fractions import Fraction as F
from experiments.marginal_projector_extendibility import replay,_telescoping_diagonal
from experiments.marginal_local_hubbard_block import _reflection


Y={102:1,612:-1}


def local_value(state):
    digits=[(state//4**i)%4 for i in range(6)]
    first=sum(digits[i]*4**i for i in range(5))
    last=sum(digits[i+1]*4**i for i in range(5))
    return Y.get(first,0)-Y.get(last,0)


class DiagonalTelescopingTests(unittest.TestCase):
    def certificate(self):
        return {'kind':'hubbard_projector_extension_v5','chain_sites':10,
                'target':{'U':'0','t':'0','V':'0'},
                'local_window':{'kind':'local_hubbard_block_v1','sites':6,'U':'0','t':'0','V':'0'},
                'vector':{0x999:1,0x666:1},'windows':2,'projector_sum_ceiling':'2',
                'penalty':'0','penalized_lower':str(min(local_value(s) for s in range(4096))),
                'joint':{'vector':{62:1,3008:1},'windows':2,'ratio':'1/2',
                         'projector_sum_ceiling':'2','penalty':'0'},
                'telescoping_diagonal':Y}

    def test_all_fock_diagonal_minimum_and_reflection(self):
        c=self.certificate();result=replay(c)
        self.assertEqual(result['local_sum_dimensions'],4096)
        self.assertEqual(result['local_maximum_psd_dimension'],200)
        self.assertEqual(F(result['periodic_lower_density']),F(c['penalized_lower'])/5)
        self.assertEqual(result['telescoping_nonzero_entries'],2)
        for s in range(4096):self.assertEqual(local_value(s),local_value(_reflection(s,6)[0]))
        c['penalized_lower']='0'
        with self.assertRaises(ValueError):replay(c)

    def test_periodic_cancellation_on_all_eight_site_determinants(self):
        # Independent cyclic base-four embedding, including windows that cross
        # the closing bond. Diagonal even operators have no CAR sign factor.
        powers=[4**i for i in range(6)]
        for state in range(4**8):
            digits=[(state//4**i)%4 for i in range(8)]
            total=0
            for offset in range(8):
                first=sum(digits[(offset+i)%8]*powers[i] for i in range(5))
                last=sum(digits[(offset+i+1)%8]*powers[i] for i in range(5))
                total+=Y.get(first,0)-Y.get(last,0)
            self.assertEqual(total,0)

    def test_refuses_malformed_or_nonodd_corrections(self):
        for source in ({102:1},{102:1,612:1},{1024:1},
                       {102:1,'102':-1},{102:1.0,612:-1.0},None,{},
                       {102:F(1,10**7),612:-F(1,10**7)}):
            with self.subTest(source=source),self.assertRaises(ValueError):_telescoping_diagonal(source)
        c=self.certificate();c['kind']='hubbard_projector_extension_v4'
        with self.assertRaisesRegex(ValueError,'requires v5'):replay(c)
        c=self.certificate();c.pop('telescoping_diagonal')
        with self.assertRaises(ValueError):replay(c)


if __name__=='__main__':unittest.main()
