"""Direct Hamiltonian MPO/DMRG proposal and charge-MPS export.

No full-state or full-H conversion is called. Numerical optimization is
number-penalty constrained; exact charge flow is enforced at export.
"""
import argparse
import importlib.metadata
import json
from pathlib import Path
import time
import numpy as np
import quimb.tensor as qtn
from quimb.operator import SparseOperatorBuilder
from experiments.marginal_symbolic import decode,encode,add,scale,mono
from research.correlated_pair_20260913.mps_round import rationalize
from research.molecular_collective_20260913.core import digest


def build_mpo(data,nalpha=None,nbeta=None,penalty=0.):
    builder=SparseOperatorBuilder(jordan_wigner=True,pauli_decompose=False,dtype=float)
    for word,c in decode(data['hamiltonian'],data['modes'],4).items():
        builder.add_term(float(c),*((('+' if creation else '-'),i) for creation,i in word))
    if penalty:
        for target,spin in ((nalpha,0),(nbeta,1)):
            sites=list(range(spin,data['modes'],2));builder.add_term(penalty*target**2)
            for a,i in enumerate(sites):
                builder.add_term(penalty*(1-2*target),('n',i))
                for j in sites[a+1:]:builder.add_term(2*penalty,('n',i),('n',j))
    return builder.build_mpo(dtype=float),len(builder._terms_raw)


def run(data,counts,bond=32,sweeps=14,penalty=4.,seed=20260913,spin=None,checkpoint=None,initial=None):
    start=time.monotonic();objective=data
    if spin is not None:
        from research.ch2_validation_20260913.certify import spin_square
        objective={**data,'hamiltonian':encode(add(decode(data['hamiltonian'],data['modes'],4),scale(spin_square(data['modes']),2),mono((),-2*spin*(spin+1))))}
    mpo,terms=build_mpo(objective,*counts,penalty)
    p=qtn.MPS_rand_state(data['modes'],bond_dim=4,phys_dim=2,dtype='float64',seed=seed)
    if initial:
        old=json.loads(Path(initial).read_text());arrays=[]
        for i,edges in enumerate(old['tensors']):
            t=np.zeros((len(old['bond_charges'][i]),2,len(old['bond_charges'][i+1])))
            for a,s,b,v in edges:t[a,s,b]=v/old['denominator']
            arrays.append(t.squeeze(axis=0) if i==0 else t.squeeze(axis=2) if i==data['modes']-1 else t)
        p=qtn.MatrixProductState(arrays,shape='lpr');p.normalize()
    dm=qtn.DMRG2(mpo,p0=p,bond_dims=[8,16,bond],cutoffs=1e-12)
    dm.opts['local_eig_ham_dense']=False
    dm.opts['local_eig_tol']=1e-5
    dm.opts['local_eig_ncv']=16
    def export_arrays(state):
        state=state.copy();state.permute_arrays('lpr');tensors=[];imax=0.
        for i in range(data['modes']):
            a=np.asarray(state[i].data);l=state.bond_size(i-1,i) if i else 1;r=state.bond_size(i,i+1) if i+1<data['modes'] else 1
            a=a.reshape(l,2,r);imax=max(imax,float(np.max(np.abs(a.imag))));tensors.append(a.real.tolist())
        if imax>1e-9:raise ValueError('Complex tensors require an explicit complex rational schema')
        return tensors,imax
    if checkpoint:
        def save_checkpoint():
            ts,_=export_arrays(dm.state);cert,stats=rationalize(data,ts,counts)
            path=Path(checkpoint);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(cert,separators=(',',':'))+'\n');tmp.replace(path)
        dm._compute_post_sweep=save_checkpoint
    converged=dm.solve(tol=1e-10,max_sweeps=sweeps,verbosity=1,sweep_sequence='RL')
    state=dm.state;physical,_=build_mpo(data)
    norm=complex(qtn.expec_TN_1D(state.H,state))
    e=complex(qtn.expec_TN_1D(state.H,physical,state))/norm
    if abs(e.imag)>1e-9:raise ValueError('Complex physical energy')
    tensors,imax=export_arrays(state)
    cert,rounding=rationalize(data,tensors,counts)
    report={'fixture_sha256':digest(data),'physical_energy_float_Ha':e.real,'penalized_energy_float_Ha':float(np.real(dm.energy)),
            'numerical_penalty_expectation_Ha':float(np.real(dm.energy))-e.real,'penalty':penalty,
            'optimizer_converged':bool(converged),'sweep_energies':list(map(lambda e:float(np.real(e)),dm.energies)),
            'requested_bond':bond,'sweeps_cap':sweeps,'seed':seed,'physical_order':list(range(data['modes'])),
            'mpo_max_bond':mpo.max_bond(),'mpo_terms':terms,'tensor_shapes':[list(np.shape(t)) for t in tensors],
            'dense_tensor_entries':sum(np.size(t) for t in tensors),'maximum_tensor_imaginary_part':imax,
            'rounding':rounding,'number_treatment':'Numerical quadratic penalties; exact charge-flow tensor export.',
            'full_state_exported':False,'full_many_body_matrix_built':False,'FCI_or_Gram_teacher_loaded':False,
            'seconds':time.monotonic()-start,'versions':{x:importlib.metadata.version(x) for x in ('quimb','numpy','scipy','cotengra','autoray')}}
    report['spin_penalty_target']=spin;report['initial_MPS_path']=initial
    return cert,report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('output');p.add_argument('--alpha',type=int,required=True);p.add_argument('--beta',type=int,required=True)
    p.add_argument('--bond',type=int,default=32);p.add_argument('--sweeps',type=int,default=14);p.add_argument('--penalty',type=float,default=4.);p.add_argument('--spin',type=int);p.add_argument('--initial');a=p.parse_args()
    data=json.loads(Path(a.fixture).read_text());out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    cert,report=run(data,[a.alpha,a.beta],a.bond,a.sweeps,a.penalty,spin=a.spin,checkpoint=str(out/'checkpoint.json'),initial=a.initial)
    (out/'state.json').write_text(json.dumps(cert,separators=(',',':'))+'\n');(out/'discovery.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
