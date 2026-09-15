import copy
from fractions import Fraction as F
import unittest

from research.side_routes_20260913.supersymmetric_chain import discover,verify


def supercharge(length):
    # Independent explicit bit-string insertion; small test sizes only.
    matrix=[[0]*(1<<length) for _ in range(1<<(length+1))]
    for code in range(1<<length):
        bits=format(code,f'0{length}b')
        for j,b in enumerate(bits):
            if b=='0':
                output=bits[:j]+'11'+bits[j+1:]
                matrix[int(output,2)][code]+=(-1)**j
    return matrix


class SupersymmetricChainTests(unittest.TestCase):
    def model(self,m=6,n=3):
        return {'modes':m,'particles':n,'hopping':['1']*(m-1),'interaction':['1']*(m-1),'fields':['0']*m}

    def test_exact_dynamic_factorization_against_full_CAR_all_sectors(self):
        from research.side_routes_20260913.supersymmetric_chain import multiply,transpose,linear
        from research.side_routes_20260913.test_positive_chain import independent
        for m in range(2,7):
            q,previous=supercharge(m),supercharge(m-1)
            actual=linear((1,multiply(transpose(q),q)),(1,multiply(previous,transpose(previous))))
            for n in range(m+1):
                model=self.model(m,n)
                states,H,*_=independent(model,{'sites':['1']*m,'bonds':['1']*(m-1)})
                # CAR reference uses site zero as the least significant bit.
                for i,si in enumerate(states):
                    row=int(format(si,f'0{m}b')[::-1],2)
                    for j,sj in enumerate(states):
                        col=int(format(sj,f'0{m}b')[::-1],2)
                        self.assertEqual(actual[row][col],H[i][j]+(m-n if i==j else 0))

    def test_large_exact_certificate_and_scaled_shifted_family(self):
        model=self.model(10000,5000)
        result=verify(model,discover(model))
        self.assertEqual(result['lower'],'-5000');self.assertEqual(result['width'],'0')
        self.assertEqual(result['many_body_states_enumerated'],0)
        model=self.model();model['hopping']=['3/2']*5;model['interaction']=['3/2']*5
        model['fields']=['2/7']*6;model['offset']='1/5'
        result=verify(model,discover(model))
        self.assertEqual(F(result['lower']),F(1,5)+F(2,7)*3-F(3,2)*3)

    def test_refuses_perturbations_false_saturation_and_tampering(self):
        model=self.model();w=discover(model)
        bad=copy.deepcopy(model);bad['interaction'][2]='10001/10000'
        with self.assertRaises(ValueError):discover(bad)
        bad=copy.deepcopy(model);bad['fields'][0]='1/10000'
        with self.assertRaises(ValueError):discover(bad)
        bad=copy.deepcopy(model);bad['extra_terms']=['unrepresented']
        with self.assertRaises(ValueError):discover(bad)
        for key,value in (('lower','0'),('upper','0'),('cycle_pair',[[0],[0],[1],[0]]),('q',[[0,0],[0,0],[0,0],[2,0]])):
            bad=copy.deepcopy(w);bad[key]=value
            with self.assertRaises(ValueError):verify(model,bad)
        model=self.model(5,2);w=discover(model);self.assertIsNone(verify(model,w)['upper'])
        w['upper']=w['lower']
        with self.assertRaises(ValueError):verify(model,w)


if __name__=='__main__':unittest.main()
