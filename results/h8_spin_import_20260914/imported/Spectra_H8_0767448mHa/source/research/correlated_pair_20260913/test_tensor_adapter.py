"""Independent tiny CAR energy oracle for Quimb MPO and paired-site adapters."""
import unittest
import numpy as np
import quimb.tensor as qtn
from research.correlated_pair_20260913.test_mps import correlated_fixture
from research.correlated_pair_20260913.mps_exact import check
from research.correlated_pair_20260913.mps_direct import build_mpo
from research.correlated_pair_20260913.mps_spatial import group_mpo,split_state
from research.correlated_pair_20260913.mps_round import rationalize


class TensorAdapter(unittest.TestCase):
    def test_grouped_mpo_and_charge_gauge(self):
        data,cert=correlated_fixture();single=[]
        for i,edges in enumerate(cert['tensors']):
            t=np.zeros((len(cert['bond_charges'][i]),2,len(cert['bond_charges'][i+1])))
            for a,s,b,v in edges:t[a,s,b]=v
            single.append(t)
        paired=[]
        for i in (0,2):
            t=np.einsum('asb,btc->astc',single[i],single[i+1]).reshape(single[i].shape[0],4,single[i+1].shape[2]);paired.append(t[0] if i==0 else t[:,:,0])
        state=qtn.MatrixProductState(paired,shape='lpr');mpo=group_mpo(build_mpo(data,1,1,4)[0])
        value=qtn.expec_TN_1D(state.H,mpo,state)/qtn.expec_TN_1D(state.H,state)
        expected=float(check(data,cert)['upper_float_Ha']);self.assertAlmostEqual(value,expected,places=12)
        rounded,stats=rationalize(data,split_state(state),[1,1]);self.assertAlmostEqual(check(data,rounded)['upper_float_Ha'],expected,places=8)


if __name__=='__main__':unittest.main()
