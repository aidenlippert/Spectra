"""Tiny independent state expansion oracles, never used by production replay."""
import copy
from fractions import Fraction as F
from itertools import product
import unittest
from research.correlated_pair_20260913.mps_exact import State,check
from research.molecular_collective_20260913.core import digest


def correlated_fixture():
    data={'modes':4,'particles':2,'hamiltonian':[
        {'word':[[1,0],[0,0]],'coefficient':'2/3'},
        {'word':[[1,0],[1,1],[0,3],[0,2]],'coefficient':'1/7'},
        {'word':[[1,2],[1,3],[0,1],[0,0]],'coefficient':'1/7'}]}
    # |0011> + 2|1100> with bit zero denoting physical mode zero.
    cert={'kind':'integer_charge_mps_v1','fixture_sha256':digest(data),'modes':4,'particles':2,
          'spin_counts':[1,1],'denominator':1,
          'bond_charges':[[[0,0]],[[0,0],[1,0]],[[0,0],[1,1]],[[1,0],[1,1]],[[1,1]]],
          'tensors':[[[0,0,0,2],[0,1,1,1]],[[0,0,0,1],[1,1,1,1]],
                     [[0,1,0,1],[1,0,1,1]],[[0,1,0,1],[1,0,0,1]]]}
    return data,cert


def explicit_moment(word):
    amps={3:1,12:2};total=0
    for state,amplitude in amps.items():
        sign=1;target=state
        for creation,i in reversed(word):
            if ((target>>i)&1)==creation: sign=0;break
            if (target&((1<<i)-1)).bit_count()%2:sign=-sign
            target^=1<<i
        total+=sign*amplitude*amps.get(target,0)
    return F(total,5)


class MPSChecks(unittest.TestCase):
    def test_all_tiny_words(self):
        data,cert=correlated_fixture();state=State(data,cert)
        self.assertEqual(state.norm_integer,5)
        letters=list(product((0,1),range(4)))
        for degree in (0,1,2,3,4):
            for word in product(letters,repeat=degree):
                self.assertEqual(state.moment(word),explicit_moment(word),word)

    def test_energy_trie(self):
        data,cert=correlated_fixture();receipt=check(data,cert)
        self.assertEqual(F(receipt['upper_Ha']),F(2,15)+F(4,35))
        self.assertEqual(receipt['stats']['enumerated_determinants'],0)

    def test_refusals(self):
        data,base=correlated_fixture()
        for change in ('charge','boundary','binding','zero','duplicate','integer','denominator'):
            cert=copy.deepcopy(base)
            if change=='charge':cert['tensors'][0][0][1]=1
            if change=='boundary':cert['bond_charges'][-1]=[[2,0]]
            if change=='binding':cert['fixture_sha256']='wrong'
            if change=='zero':cert['tensors'][1]=[]
            if change=='duplicate':cert['tensors'][0].append(cert['tensors'][0][0])
            if change=='integer':cert['tensors'][0][0][-1]=0.5
            if change=='denominator':cert['denominator']=-1
            with self.subTest(change=change),self.assertRaises(ValueError):check(data,cert)


if __name__=='__main__':unittest.main()
