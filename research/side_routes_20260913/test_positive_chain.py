import copy
import unittest
from fractions import Fraction as F
from itertools import combinations
from research.side_routes_20260913.positive_chain import bounds,local_energy,replay


def independent(model,amp):
    m,n=model['modes'],model['particles'];t=list(map(F,model['hopping']));v=list(map(F,model['interaction']));h=list(map(F,model['fields']))
    states=[]
    for indices in combinations(range(m),n):states.append(sum(1<<i for i in indices))
    pos={s:i for i,s in enumerate(states)}
    H=[[F(0) for _ in states] for _ in states];psi=[]
    for col,s in enumerate(states):
        bits=[(s>>i)&1 for i in range(m)]
        p=F(1)
        for i,b in enumerate(bits):p*=F(amp['sites'][i])**b
        for i in range(m-1):p*=F(amp['bonds'][i])**(bits[i]*bits[i+1])
        psi.append(p)
        H[col][col]=F(model.get('offset',0))+sum(h[i]*bits[i] for i in range(m))+sum(v[i]*bits[i]*bits[i+1] for i in range(m-1))
        # Independent CAR action for the two ordered nearest-neighbor hops.
        for i in range(m-1):
            for create,annihilate in ((i,i+1),(i+1,i)):
                state=s;sign=1
                for flag,j in ((0,annihilate),(1,create)):
                    if bool(state&(1<<j))==bool(flag):sign=0;break
                    if (state&((1<<j)-1)).bit_count()%2:sign=-sign
                    state^=1<<j
                if sign:H[pos[state]][col]-=t[i]*sign
    energies=[sum(H[i][j]*psi[j] for j in range(len(states)))/psi[i] for i in range(len(states))]
    Z=sum(p*p for p in psi);upper=sum(p*p*e for p,e in zip(psi,energies))/Z
    return states,H,psi,energies,Z,upper


class PositiveChainTests(unittest.TestCase):
    def test_exact_DP_against_independent_CAR_every_sector(self):
        m=6;model={'modes':m,'particles':0,'hopping':['1','2/3','3/5','5/7','7/11'],
            'interaction':['1/3','-2/5','2/7','3/4','-1/2'],'fields':['1/7','-1/3','2/5','-3/7','4/9','-5/11'],'offset':'2/13'}
        amp={'sites':['1','3/2','4/3','5/4','6/5','7/6'],'bonds':['2/3','3/4','4/5','5/6','6/7']}
        for n in range(m+1):
            model['particles']=n;states,H,psi,energies,Z,upper=independent(model,amp);r=bounds(model,amp)
            self.assertEqual(F(r['lower']),min(energies));self.assertEqual(F(r['upper']),upper);self.assertEqual(F(r['partition']),Z)
            self.assertEqual(local_energy(model,amp,list(map(int,r['minimum_configuration']))),min(energies))
            for state,e in zip(states,energies):self.assertEqual(local_energy(model,amp,[(state>>i)&1 for i in range(m)]),e)

    def test_nonfree_frustration_free_control_large_chain(self):
        m=64;model={'modes':m,'particles':m//2,'hopping':['1']*(m-1),'interaction':['-2']*(m-1),'fields':['1']+['2']*(m-2)+['1']}
        amp={'sites':['1']*m,'bonds':['1']*(m-1)};r=bounds(model,amp)
        self.assertEqual(F(r['lower']),0);self.assertEqual(F(r['upper']),0)
        self.assertLessEqual(r['peak_dp_states'],8*(m//2+1));self.assertEqual(r['many_body_states_enumerated'],0)

    def test_numerical_spectrum_between_exact_bounds(self):
        import numpy as np
        m=6;model={'modes':m,'particles':3,'hopping':['1']*(m-1),'interaction':['1/2']*(m-1),'fields':['0','1/3','-1/5','2/7','0','-1/4']}
        amp={'sites':['1','4/5','6/5','1','1','5/4'],'bonds':['3/4']*(m-1)}
        _,H,*_=independent(model,amp);energy=np.linalg.eigvalsh(np.array(H,dtype=float))[0];r=bounds(model,amp)
        self.assertLessEqual(r['lower_float'],energy);self.assertGreaterEqual(r['upper_float'],energy)

    def test_uniform_parent_perturbation_has_extensive_width(self):
        # Independent combinatorics: E_local=delta times occupied adjacent pairs.
        # Their uniform-sector expectation is N(N-1)/M; minimum is zero at half filling.
        delta=F(1,10000)
        for m in (4,8,16,32,64):
            n=m//2
            model={'modes':m,'particles':n,'hopping':['1']*(m-1),
                   'interaction':[str(-2+delta)]*(m-1),'fields':['1']+['2']*(m-2)+['1']}
            amp={'sites':['1']*m,'bonds':['1']*(m-1)}
            r=bounds(model,amp)
            self.assertEqual(F(r['lower']),0)
            self.assertEqual(F(r['width']),delta*F(n*(n-1),m))

    def test_refusal_and_claim_tampering(self):
        model={'modes':2,'particles':1,'hopping':['1'],'interaction':['0'],'fields':['0','0']};amp={'sites':['1','1'],'bonds':['1']}
        r=bounds(model,amp);payload={'model':model,'amplitude':amp,'claim':r};replay(payload)
        bad=copy.deepcopy(payload);bad['claim']['upper']='-999'
        with self.assertRaises(ValueError):replay(bad)
        bad=copy.deepcopy(model);bad['hopping']=['-1']
        with self.assertRaises(ValueError):bounds(bad,amp)
        bad=copy.deepcopy(amp);bad['sites']=['0','1']
        with self.assertRaises(ValueError):bounds(model,bad)
        bad=copy.deepcopy(model);bad['particles']=True
        with self.assertRaises(ValueError):bounds(bad,amp)
        bad=copy.deepcopy(amp);bad['sites']=[1.0,1]
        with self.assertRaises(ValueError):bounds(model,bad)

if __name__=='__main__':unittest.main()
