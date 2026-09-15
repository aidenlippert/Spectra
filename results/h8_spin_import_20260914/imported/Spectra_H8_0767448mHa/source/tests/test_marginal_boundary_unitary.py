import json
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_boundary_unitary import (
    rotation,energy_shift,edge_reductions,shift_polynomial,family_bound,replay,model_shift,replay_model)
from experiments.marginal_hopping_filter import boundary_data
from tests.test_marginal_hopping_filter import VALENCE,action,hopping,dot

ROOT=Path(__file__).resolve().parents[1]


def reduced(vector,offset):
    norm=dot(vector,vector)
    groups={}
    for s,a in vector.items():
        environment=(s&((1<<offset)-1))|((s>>(offset+4))<<offset)
        groups.setdefault(environment,{})[(s>>offset)&15]=a
    rho=[[F(0)]*16 for _ in range(16)]
    for part in groups.values():
        for i,a in part.items():
            for j,b in part.items():rho[i][j]+=F(a*b,norm)
    return rho


class BoundaryUnitaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs=json.loads((ROOT/'results/marginal_graded_hubbard8/tiled_eight_upper/certificate.json').read_text())
        cls.data=edge_reductions(cls.inputs['upper'],cls.inputs['hamiltonian'])

    def test_rotation_and_identity(self):
        self.assertEqual(rotation(0),{i:{i:F(1)} for i in range(16)})
        gate=rotation('1/10')
        for j,col in gate.items():self.assertEqual(sum(a*a for a in col.values()),1)
        self.assertEqual(energy_shift(self.data,0),0)
        with self.assertRaises(ValueError):rotation(.1)

    def test_two_and_three_blocks_direct_full_hubbard_contraction(self):
        data={'left':reduced(VALENCE,0),'right':reduced(VALENCE,4)}
        shift=energy_shift(data,'1/10')
        gate=rotation('1/10')
        # 101 times the unitary is an integer matrix; use exact integer states.
        gate={i:{j:int(a*101) for j,a in col.items()} for i,col in gate.items()}
        for q in (2,3):
            with self.subTest(blocks=q):
                sites=4*q;state={0:1}
                for block in range(q):
                    state={s+(t<<(8*block)):a*b for s,a in state.items() for t,b in VALENCE.items()}
                for cut in range(4,sites,4):
                    offset=2*(cut-1);mask=15<<offset;out={}
                    for s,a in state.items():
                        for local,b in gate[(s>>offset)&15].items():
                            target=(s&~mask)|(local<<offset);out[target]=out.get(target,0)+a*b
                    state={s:a for s,a in out.items() if a}
                    self.assertLessEqual(len(state),65536)
                norm=dot(state,state)
                self.assertEqual(norm,12**q*101**(2*(q-1)))
                onsite=sum(a*a*4*sum(((s>>(2*i))&3)==3 for i in range(sites)) for s,a in state.items())
                kinetic=sum(dot(state,action(hopping(cut),state)) for cut in range(1,sites))
                self.assertEqual(F(onsite+kinetic,norm),(q-1)*shift)
                self.assertEqual(reduced(state,0),data['left'])
                self.assertEqual(reduced(state,2*(sites-2)),data['right'])
                for edge in (0,sites-1):
                    self.assertEqual(sum(a*a*int(((s>>(2*edge))&3)==3) for s,a in state.items()),0)

    def test_h8_reduced_observables_match_independent_stream(self):
        old=boundary_data(self.inputs['upper'],self.inputs['hamiltonian'])
        for side,site,index in [('left',0,0),('right',1,1)]:
            rho=[[F(x) for x in row] for row in self.data[side]]
            self.assertEqual(sum(rho[i][i] for i in range(16)),1)
            doublons=sum(rho[i][i] for i in range(16) if ((i>>(2*site))&3)==3)
            self.assertEqual(doublons,F(old['edge_doublons'][index]))
            value=F(0)
            for s in range(16):
                for t,a in action(hopping(1,1),{s:1}).items():value+=a*rho[s][t]
            self.assertEqual(value,F(old['edge_hopping'][index]))

    def test_global_quartic_certificate_and_rejects_stronger_claim(self):
        bound=family_bound(self.data,'-0.110426865','-0.13881')
        achieved=energy_shift(self.data,'0.150424')
        self.assertTrue(0<=achieved-F(bound['minimum_shift_lower'])<F(1,10**9))
        coefficients=list(map(F,bound['numerator_coefficients']))
        for x in (F(1,3),F(-7,5)):
            self.assertEqual(sum(a*x**k for k,a in enumerate(coefficients))/(1+x*x)**2,
                             energy_shift(self.data,x))
        with self.assertRaises(ValueError):family_bound(self.data,'-0.110426','-0.13881')

    def test_physical_replay_and_fixed_size_arithmetic(self):
        result=replay(self.inputs['upper'],self.inputs['hamiltonian'],'0.150424',1000000)
        self.assertEqual(result['norm_factor'],'1')
        self.assertEqual(F(result['upper']),125000*F(self.data['energy8'])+124999*energy_shift(self.data,'0.150424'))
        with self.assertRaises(ValueError):replay(self.inputs['upper'],self.inputs['hamiltonian'],'0.1',7)
        with self.assertRaises(ValueError):replay(self.inputs['upper'],self.inputs['hamiltonian'],'0.1',18)

    def test_transfer_to_density_interaction_by_direct_three_block_state(self):
        data={'left':reduced(VALENCE,0),'right':reduced(VALENCE,4)}
        q,sites=3,12
        for operation in ('unitary','linear_filter'):
            with self.subTest(operation=operation):
                if operation=='unitary':
                    gate={i:{j:int(a*101) for j,a in col.items()} for i,col in rotation('1/10').items()}
                    parameter,denominator='1/10',101
                else:
                    gate={}
                    for i in range(16):
                        h=action(hopping(1),{i:1})
                        gate[i]={j:5*int(i==j)-h.get(j,0) for j in range(16)
                                 if 5*int(i==j)-h.get(j,0)}
                    parameter,denominator='1/5',5
                state={0:1}
                for block in range(q):
                    state={s+(t<<(8*block)):a*b for s,a in state.items() for t,b in VALENCE.items()}
                for cut in (4,8):
                    offset=2*(cut-1);mask=15<<offset;out={}
                    for s,a in state.items():
                        for local,b in gate[(s>>offset)&15].items():
                            target=(s&~mask)|(local<<offset);out[target]=out.get(target,0)+a*b
                    state={s:a for s,a in out.items() if a}
                norm=dot(state,state)
                self.assertEqual(norm,12**q*(101**2 if operation=='unitary' else 26)**(q-1))
                onsite=sum(a*a*3*sum(((s>>(2*i))&3)==3 for i in range(sites)) for s,a in state.items())
                kinetic=F(2,3)*sum(dot(state,action(hopping(cut),state)) for cut in range(1,sites))
                density=F(0)
                for s,a in state.items():
                    charge=[((s>>(2*i))&3).bit_count()-1 for i in range(sites)]
                    density+=F(a*a,2)*sum(charge[i]*charge[i+1] for i in range(sites-1))
                self.assertNotEqual(density,0)
                shift=model_shift(data,parameter,operation,U='3',t='2/3',V='1/2')
                self.assertEqual((onsite+kinetic+density)/norm,(q-1)*shift)
                self.assertEqual(reduced(state,0),data['left'])
                self.assertEqual(reduced(state,20),data['right'])

    def test_transferred_source_energy_and_target_refusals(self):
        reference=model_shift(self.data,'57277/250000','linear_filter')
        old=boundary_data(self.inputs['upper'],self.inputs['hamiltonian'])
        from experiments.marginal_hopping_filter import merge_shift
        self.assertEqual(reference,merge_shift(old['delta'],'57277/250000'))
        r=replay_model(self.inputs['upper'],self.inputs['hamiltonian'],'57277/250000',24,
                       operation='linear_filter',U='3',t='2/3',V='1/2')
        self.assertEqual(r['target'],{'U':'3','t':'2/3','V':'1/2'})
        D=F(self.data['total_doublons']);T=4*D-F(self.data['energy8'])
        self.assertEqual(F(r['target_block_energy']),3*D-F(2,3)*T+F(self.data['internal_density_correlation'])/2)
        for kwargs in ({'sites':18},{'sites':24,'V':.5},{'sites':24,'operation':'bad'}):
            with self.assertRaises(ValueError):
                replay_model(self.inputs['upper'],self.inputs['hamiltonian'],'1/10',**kwargs)


if __name__=='__main__':unittest.main()
