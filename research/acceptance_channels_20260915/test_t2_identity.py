from fractions import Fraction
from itertools import combinations
import unittest
from experiments.marginal_symbolic import add, adj, canonical, mono, scale
from research.acceptance_channels_20260915.select_channels import t2_terms
from research.acceptance_channels_20260915.augment import paired_entry


class T2IdentityTest(unittest.TestCase):
    def test_all_three_orbital_cross_terms_against_literal_CAR(self):
        triples=[(p,q,r) for p in range(3) for q,r in combinations(range(3),2)]
        for p,q,r in triples:
            for s,t,u in triples:
                raw=add(*(mono(w,c) for w,c in t2_terms(p,q,r,s,t,u)))
                actual=canonical(scale(add(raw,adj(raw)),Fraction(1,2)))
                expected=paired_entry(((1,2*p),(0,2*q+1),(0,2*r+1)),((1,2*s),(0,2*t+1),(0,2*u+1)))
                self.assertEqual(actual,expected)


if __name__=='__main__':unittest.main()
