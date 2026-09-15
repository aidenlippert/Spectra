"""The same direct MPS proposal with paired spin sites in each DMRG update."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import quimb.tensor as qtn
from research.correlated_pair_20260913.mps_direct import build_mpo
from research.correlated_pair_20260913.mps_round import rationalize
from research.molecular_collective_20260913.core import digest


def group_mpo(mpo):
    mpo=mpo.copy();mpo.permute_arrays('lrud');arrays=[]
    for i in range(0,mpo.L,2):
        a=np.asarray(mpo[i].data).reshape(mpo.bond_size(i-1,i) if i else 1,mpo.bond_size(i,i+1),2,2)
        b=np.asarray(mpo[i+1].data).reshape(mpo.bond_size(i,i+1),mpo.bond_size(i+1,i+2) if i+2<mpo.L else 1,2,2)
        t=np.einsum('abuv,bcxy->acuxvy',a,b).reshape(a.shape[0],b.shape[1],4,4)
        arrays.append(t[0] if i==0 else t[:,0] if i+2==mpo.L else t)
    return qtn.MatrixProductOperator(arrays,shape='lrud')


def split_state(state):
    state=state.copy();state.permute_arrays('lpr');tensors=[]
    for i in range(state.L):
        l=state.bond_size(i-1,i) if i else 1;r=state.bond_size(i,i+1) if i+1<state.L else 1
        a=np.asarray(state[i].data).reshape(l*2,2*r)
        u,s,v=np.linalg.svd(a,full_matrices=False);rank=len(s)
        tensors.extend([u.reshape(l,2,rank).tolist(),(s[:,None]*v).reshape(rank,2,r).tolist()])
    return tensors


def run(data,counts,out,bond=64,sweeps=10,initial=None):
    start=time.monotonic();raw,terms=build_mpo(data,*counts,4.);mpo=group_mpo(raw)
    before_compression=mpo.max_bond();mpo.compress(cutoff=1e-12)
    # Two-spatial-orbital local updates can create pair correlations directly.
    # This exact HF product seed has the requested Nalpha,Nbeta.
    bits=[(2 if i<counts[0] else 0)+(1 if i<counts[1] else 0) for i in range(data['modes']//2)]
    arrays=[np.eye(4)[v].reshape((1,4,1)) for v in bits];arrays[0]=arrays[0][0];arrays[-1]=arrays[-1][:,:,0]
    seed=qtn.MPS_rand_state(data['modes']//2,bond_dim=4,phys_dim=4,dtype='float64',seed=20260913)
    if initial:
        old=json.loads(Path(initial).read_text());single=[]
        if old['fixture_sha256']!=digest(data):raise ValueError('Initial MPS fixture mismatch')
        for i,edges in enumerate(old['tensors']):
            t=np.zeros((len(old['bond_charges'][i]),2,len(old['bond_charges'][i+1])))
            for a,s,b,v in edges:t[a,s,b]=v/old['denominator']
            single.append(t)
        arrays=[]
        for i in range(0,data['modes'],2):
            t=np.einsum('asb,btc->astc',single[i],single[i+1]).reshape(single[i].shape[0],4,single[i+1].shape[2])
            arrays.append(t[0] if i==0 else t[:,:,0] if i+2==data['modes'] else t)
        seed=qtn.MatrixProductState(arrays,shape='lpr');seed.normalize()
    dm=qtn.DMRG2(mpo,p0=seed,bond_dims=[bond] if initial else [8,16,bond],cutoffs=1e-12);dm.opts['local_eig_ham_dense']=False;dm.opts['local_eig_tol']=1e-5;dm.opts['local_eig_ncv']=16
    def checkpoint():
        cert,stats=rationalize(data,split_state(dm.state),counts)
        (out/'checkpoint.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    dm._compute_post_sweep=checkpoint
    converged=dm.solve(tol=1e-9,max_sweeps=sweeps,sweep_sequence='RL',verbosity=1)
    state=dm.state;physical=group_mpo(build_mpo(data)[0]);energy=float(qtn.expec_TN_1D(state.H,physical,state)/qtn.expec_TN_1D(state.H,state))
    tensors=split_state(state);cert,stats=rationalize(data,tensors,counts)
    report={'fixture_sha256':digest(data),'physical_energy_float_Ha':energy,'penalized_energy_float_Ha':float(dm.energy),
            'initial_MPS_path':initial,
            'requested_spatial_bond':bond,'spatial_bonds':[state.bond_size(i,i+1) for i in range(state.L-1)],
            'export_spin_bonds':[len(q) for q in cert['bond_charges']],'export_nonzero_entries':sum(map(len,cert['tensors'])),
            'MPO_bond':mpo.max_bond(),'MPO_bond_before_compression':before_compression,'MPO_compression_cutoff_proposal_only':1e-12,'sweeps_cap':sweeps,'sweep_energies':list(map(float,dm.energies)),
            'converged':bool(converged),'rounding':stats,'seconds':time.monotonic()-start,
            'initializer':'Random MPS bond4, seed20260913, paired spin sites; no FCI or Gram teacher','full_state_or_H_expansion':False}
    (out/'state.json').write_text(json.dumps(cert,separators=(',',':'))+'\n');(out/'discovery.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('output',type=Path);p.add_argument('--alpha',type=int,required=True);p.add_argument('--beta',type=int,required=True);p.add_argument('--bond',type=int,default=64);p.add_argument('--sweeps',type=int,default=10);p.add_argument('--initial');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);run(json.loads(Path(a.fixture).read_text()),[a.alpha,a.beta],a.output,a.bond,a.sweeps,a.initial)
