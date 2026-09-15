"""Hamiltonian-only CISD witness generation and exact rational replay."""
from pathlib import Path
import json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from experiments.marginal_determinant_tree import DeterminantOracle
from research.certificate_scaling.streaming_reference_upper import upper as stream_upper

def run(fixture, out, scale=10**12):
    t=time.monotonic(); cert=json.loads(fixture.read_text()); oracle=DeterminantOracle(cert)
    states=[s for s in range(1<<oracle.modes) if oracle.valid_state(s)]
    # Reference selected solely from diagonal H; no saved wavefunction is read.
    def diagonal_energy(s):
        return sum(float(c) for support,c in oracle.diagonal.items() if s & support == support)
    diag=[(diagonal_energy(s),s) for s in states]
    hf=min(diag)[1]; occ=[i for i in range(oracle.modes) if hf>>i&1]; vir=[i for i in range(oracle.modes) if not hf>>i&1]
    oracle.cache.clear()  # diagonal scan is discovery only; replay has its own bounded cache
    basis={hf};
    for i in occ:
        for a in vir: basis.add(hf^(1<<i)^(1<<a))
    for x,i in enumerate(occ):
        for j in occ[x+1:]:
            for y,a in enumerate(vir):
                for b in vir[y+1:]: basis.add(hf^(1<<i)^(1<<j)^(1<<a)^(1<<b))
    basis=sorted(basis); n=len(basis); mat=np.zeros((n,n))
    for j,s in enumerate(basis):
        for target,v in oracle.action(s).items():
            if target in basis: mat[basis.index(target),j]=float(v)
    ev,vec=np.linalg.eigh(mat); coeff=vec[:,0]; amps=[int(round(float(x)*scale)) for x in coeff]
    witness={'states':[s for s,a in zip(basis,amps) if a], 'amplitudes':[a for a in amps if a]}
    val,stats=stream_upper(cert,witness)
    out.mkdir(parents=True,exist_ok=True); (out/'upper.json').write_text(json.dumps({'independent_upper':witness,'upper':str(val)},indent=2)+'\n')
    rec={'kind':'hamiltonian_only_cisd_witness_v1','fixture':str(fixture),'reference_state':hf,'cisd_basis_size':n,'witness_support':len(witness['states']),'amplitude_scale':scale,'upper':str(val),'upper_float':float(val),'discovery_seconds':time.monotonic()-t,'matrix_entries':n*n,'upper_replay':stats,'discovery_inputs':['fixture Hamiltonian only'],'validation_note':'FCI/scoreboard values are not used for discovery; compare externally after freeze. CISD means at most two particle-hole substitutions from lowest-diagonal determinant.'}
    (out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n'); return rec

if __name__=='__main__':
    root=ROOT/'results/certificate_scaling/active_space_ladder'; outroot=ROOT/'results/all_angles_20260913/upper_states'; rows=[]
    for n in (4,6,8):
        try: rows.append(run(root/f'h{n}/fixture.json',outroot/f'h{n}'))
        except Exception as e: rows.append({'system':f'H{n}','error':repr(e)})
    (outroot/'summary.json').write_text(json.dumps(rows,indent=2)+'\n'); print(json.dumps(rows,indent=2))
