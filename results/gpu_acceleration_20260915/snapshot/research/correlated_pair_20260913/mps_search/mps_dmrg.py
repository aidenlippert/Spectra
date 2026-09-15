"""Direct CAR-to-MPO DMRG proposer using isolated Quimb.

This path deliberately never constructs a Fock-space matrix or FCI vector.
Particle number is currently enforced by the initial product state only; the
receipt labels the run penalty-constrained until an explicit U(1) HilbertSpace
sector is added.
"""
from pathlib import Path
raise RuntimeError('Superseded incomplete tensor proposer: use mps_direct or mps_spatial with mps_exact acceptance.')
import json, sys, time
import numpy as np
from quimb.operator import SparseOperatorBuilder
import quimb.tensor as qtn
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_symbolic import decode

ROOT=Path(__file__).resolve().parents[3]

def build_mpo(data, nalpha=None, nbeta=None, penalty=0.0):
    m=data['modes']; b=SparseOperatorBuilder(jordan_wigner=True, pauli_decompose=True, dtype=complex)
    for w,c in decode(data['hamiltonian'],m,4).items():
        b.add_term(float(c), *((('+' if int(op) else '-'), int(i)) for op,i in w))
    if penalty:
        for target, spin in ((nalpha,0),(nbeta,1)):
            sites=list(range(spin,m,2))
            for i in sites: b.add_term(penalty*(1-2*target), ('n',i))
            for a,i in enumerate(sites):
                for j in sites[a+1:]: b.add_term(2*penalty, ('n',i), ('n',j))
            b.add_term(penalty*target*target)
    return b.build_mpo(dtype=complex), len(b._terms_raw)

def product_state(m, nalpha, nbeta):
    # Interleaved spin orbital order; this is only an HF/product initializer.
    bits=[1 if (i//2 < nalpha and i%2==0) or (i//2 < nbeta and i%2==1) else 0 for i in range(m)]
    return qtn.MPS_computational_state(''.join(map(str,bits)))

def optimize(data, nalpha, nbeta, bond=16, max_sweeps=8, seed=20260913, penalty=20.0):
    start=time.monotonic(); mpo,nterms=build_mpo(data,nalpha,nbeta,penalty); rng=np.random.default_rng(seed)
    p=qtn.MPS_rand_state(data['modes'], bond_dim=min(4,bond), phys_dim=2, dtype='float64', seed=seed)
    dm=qtn.DMRG2(mpo, p0=p, bond_dims=[min(8,bond), bond], cutoffs=1e-8)
    dm.solve(max_sweeps=max_sweeps, verbosity=0)
    mps=dm.state
    physical_mpo,_=build_mpo(data)
    physical_energy=float(np.real(mps.expec(physical_mpo).data.ravel()[0]))
    def tj(t):
        a=np.asarray(t.data).ravel(); return [float(z.real) if abs(z.imag)<1e-10 else [float(z.real),float(z.imag)] for z in a]
    return {'energy_penalized_float':float(np.real(dm.energy)), 'energy_float_unpenalized':physical_energy, 'bond_dimension':bond, 'sweeps':max_sweeps,
            'tensor_shapes':[list(t.shape) for t in mps.tensors], 'tensor_entries_flat':[tj(t) for t in mps.tensors], 'physical_mode_order':list(range(data['modes'])), 'hamiltonian_terms':nterms,
            'seconds':time.monotonic()-start, 'particle_sector':'HF initializer only; U(1) not yet certified',
            'full_fock_matrix_built':False,'full_state_exported':False,'seed':seed,'penalty':penalty}, mps

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('fixture'); ap.add_argument('--na',type=int,required=True); ap.add_argument('--nb',type=int,required=True); ap.add_argument('--bond',type=int,default=16)
    a=ap.parse_args(); d=json.loads(Path(a.fixture).read_text()); row,_=optimize(d,a.na,a.nb,a.bond); print(json.dumps(row,indent=2))
