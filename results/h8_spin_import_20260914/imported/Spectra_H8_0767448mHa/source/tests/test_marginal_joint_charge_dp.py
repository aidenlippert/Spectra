import copy
from fractions import Fraction as F
from itertools import combinations
import unittest

from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_joint_charge_dp import joint_envelope,replay


def example():
    h=build(4,4,1)
    return dict(kind='joint_charge_product_dp_v1',sites=4,U='4',t='1',target_lower='-100',
                **{k:h[k] for k in ('modes','particles','hamiltonian')},
                onsite=[['1','2','3'],['2','1','1'],['1','3','2'],['3','2','1']],
                pairs=[{'i':0,'j':1,'values':[['1','2','1'],['3','1','2'],['2','3','1']]},
                       {'i':1,'j':3,'values':[['2','1','3'],['1','2','1'],['3','1','2']]}])


class JointChargeDPTests(unittest.TestCase):
    def test_joint_maximum_equals_independent_fermionic_rows(self):
        self.compare_rows(example())

    def test_general_count_profile_matches_independent_fermionic_rows(self):
        c=example();c['doublon_weights']=['0','3/2','2/3']
        self.compare_rows(c)
        for values in [['1','1','2'],['0','1','0'],['0','1']]:
            c['doublon_weights']=values
            with self.assertRaises(ValueError):replay(c)

    def compare_rows(self,c):
        oracle=DeterminantOracle(c);floors={};pattern_floors={}
        def metric(state):
            q=[((state>>(2*i))&1)+((state>>(2*i+1))&1)-1 for i in range(4)]
            D=sum(x*x for x in q)//2;v=F(c.get('doublon_weights',['0','1','2'])[D])
            for i,x in enumerate(q):v*=F(c['onsite'][i][x+1])
            for pair in c['pairs']:v*=F(pair['values'][q[pair['i']]+1][q[pair['j']]+1])
            return D,v
        for up in combinations(range(4),2):
            for down in combinations(range(4),2):
                state=sum(1<<(2*i) for i in up)+sum(1<<(2*i+1) for i in down)
                D,v=metric(state)
                if not D:continue
                action=oracle.action(state)
                lower=action.get(state,F(0))-sum(abs(value)*metric(target)[1]/v for target,value in action.items() if target!=state)
                floors[D]=min(floors.get(D,lower),lower)
                pattern=tuple(((state>>(2*i))&1)+((state>>(2*i+1))&1)-1 for i in range(4))
                pattern_floors[pattern]=min(pattern_floors.get(pattern,lower),lower)
        result=joint_envelope(c,return_witnesses=True)
        self.assertEqual({row['doublons']:F(row['row_lower']) for row in result['rows']},floors)
        for row in result['rows']:
            q=row['worst_charge_pattern']
            self.assertEqual(sum(q),0)
            self.assertEqual(sum(x*x for x in q),2*row['doublons'])
            self.assertEqual(pattern_floors[tuple(q)],F(row['row_lower']))
        c['target_lower']=result['minimum_lower']
        self.assertTrue(replay(c)['exact_accepted'])
        c['target_lower']=str(F(result['minimum_lower'])+F(1,1000))
        with self.assertRaisesRegex(ValueError,'misses'):replay(c)

    def test_uniform_metric_and_refusal_gates(self):
        c=example();c['onsite']=[['1']*3 for _ in range(4)];c['pairs']=[]
        result=joint_envelope(c)
        self.assertEqual(result['range'],0)
        self.assertLess(result['largest_local_window'],4)
        c['target_lower']=result['minimum_lower'];self.assertTrue(replay(c)['exact_accepted'])
        for mutate in [lambda v:v['onsite'][0].__setitem__(0,'0'),
                       lambda v:v['hamiltonian'][0].__setitem__('coefficient','123'),
                       lambda v:v.__setitem__('particles',3),
                       lambda v:v.__setitem__('sites',True)]:
            bad=copy.deepcopy(c);mutate(bad)
            with self.assertRaises(ValueError):replay(bad)
        with self.assertRaisesRegex(ValueError,'state budget'):joint_envelope(c,max_states=1)

    def test_eight_site_uniform_fixed_metric_limit(self):
        h=build(8,4,1)
        c=dict(kind='joint_charge_product_dp_v1',sites=8,U='4',t='1',target_lower='-18',
               **{k:h[k] for k in ('modes','particles','hamiltonian')},onsite=[['1']*3 for _ in range(8)],pairs=[])
        receipt=replay(c)
        self.assertEqual([F(row['row_lower']) for row in receipt['rows']],[F(-18),F(-5),F(5,3),F(11,2)])


if __name__=='__main__':unittest.main()
