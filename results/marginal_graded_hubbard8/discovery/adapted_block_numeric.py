"""Bounded floating-point proposal; exact replay is the acceptance gate."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_boundary_unitary import _physical_source
from experiments.marginal_adapted_block import effective_model
from experiments.marginal_symmetry_moments import SymmetryMomentOracle

def main():
    c=json.loads((ROOT/'results/marginal_graded_hubbard8/tiled_eight_upper/certificate.json').read_text())
    original,source,_=_physical_source(c['upper'],c['hamiltonian'])
    reps=sorted({e[0] for s in range(65536) if original.valid_state(s) and (e:=original.orbit(s)) is not None})
    index={s:i for i,s in enumerate(reps)}
    weights=np.sqrt([original.orbit(s)[2] for s in reps])
    def matrix(oracle):
        rows,cols,values=[],[],[]
        for j,s in enumerate(reps):
            for t,a in oracle.action(s).items():
                i=index[t];rows.append(i);cols.append(j);values.append(float(a)*weights[i]/weights[j])
        return csc_matrix((values,(rows,cols)),shape=(len(reps),len(reps)))
    H=matrix(original)
    half=matrix(SymmetryMomentOracle(effective_model(1)[0]));B=2*(half-H)
    v=np.array([float(source.get(s,0)) for s in reps])*weights;v/=np.linalg.norm(v)
    def proposal(eta,degree=8):
        w=eta*eta/(1+eta*eta);K=H+w*B
        basis=[v];polys=[np.eye(1,degree+1,0)[0]]
        for j in range(degree):
            z=K@basis[-1]+4*basis[-1]
            p=np.roll(polys[-1],1);p[0]=0
            for _ in range(2):
                for u,up in zip(basis,polys):
                    a=u@z;z-=a*u;p-=a*up
            norm=np.linalg.norm(z)
            if norm<1e-12:break
            basis.append(z/norm);polys.append(p/norm)
        Q=np.array(basis).T
        vals,vecs=np.linalg.eigh(Q.T@(K@Q))
        p=vecs[:,0]@np.array(polys);p/=np.max(np.abs(p))
        return (vals[0]+4*w-2*eta/(1+eta*eta))/8,p
    opt=minimize_scalar(lambda eta:proposal(eta)[0],bounds=(.1,.4),method='bounded',options={'xatol':1e-10,'maxiter':40})
    eta=round(float(opt.x),6);value,p=proposal(eta)
    out={'eta':str(eta),'coefficients':[format(x,'.12f') for x in p],
         'numeric_objective_density':float(value),'orbit_matrix_dimension':len(reps),
         'scope':'Floating-point finite H8 Krylov proposal only; no accepted energy certificate.'}
    dest=ROOT/'results/marginal_graded_hubbard8/adapted_block';dest.mkdir(exist_ok=True)
    (dest/'numeric_proposal.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)

if __name__=='__main__':main()
