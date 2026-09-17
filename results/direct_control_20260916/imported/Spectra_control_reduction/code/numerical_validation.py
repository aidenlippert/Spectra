"""Independent numerical dynamics checks; not the accepting proof."""
from propose import *
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm
import random
start=time.monotonic();labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');Di=sparse.diags(D);rng=np.random.default_rng(20164)
records=[]
for name in ['short24','long64','baseline32','short8_refusal']:
 case=json.loads((ROOT/'inputs'/f'{name}.json').read_text());J=np.array(case['basis'],float)/case['basis_denominator'];zd=case['trajectory_denominator']
 actual=psi.astype(complex);dt=.25;prev=None;groups=[]
 for seg in case['segments']:
  key=(float(F(seg['u'])),float(F(seg['v'])))
  if groups and groups[-1][0]==key:groups[-1][1]+=dt
  else:groups.append([key,dt])
 for (u,v),t in groups:
  A=K+u*Di+v*W;actual=expm_multiply(-1j*t*A,actual,traceA=-1j*t*A.diagonal().sum())
 c=case['segments'][-1];z=(np.sum(c['re'],axis=0)+1j*np.sum(c['im'],axis=0))/zd;pred=J@z
 records.append({'case':name,'full_numerical_D':float(np.real(actual.conj()@(D*actual))),'full_numerical_state_norm':float(np.linalg.norm(actual)),
   'numerical_state_error':float(np.linalg.norm(pred-actual)),'basis_dimension':J.shape[1]})
 print(records[-1],flush=True)
# Prescribed corners of the long-control waveform band. No control refitting.
for du,dv in [(-.00009,-.00009),(-.00009,.00009),(.00009,-.00009),(.00009,.00009)]:
 A=K+(.2+du)*Di+(.1+dv)*W;actual=expm_multiply(-21j*A,psi,traceA=-21j*A.diagonal().sum())
 records.append({'case':'long64_perturbed','du':du,'dv':dv,'full_numerical_D':float(np.real(actual.conj()@(D*actual)))})
# 2-level noisy oracle independently checks the norm-contraction robustness bound.
X=np.array([[0.,1.],[1.,0.]]);Z=np.diag([1.,-1.]);Id=np.eye(2);rho=np.diag([1.,0.]);outnoise=[]
for gamma in [0,.001,.01,.1]:
 A=.4*X+.2*Z
 L=-1j*(np.kron(Id,A)-np.kron(A.T,Id))+gamma*(np.kron(Z,Z)-np.eye(4))
 r=expm(3*L)@rho.reshape(-1,order='F');noisy=r.reshape(2,2,order='F')
 U=expm(-3j*A);clean=U@rho@U.conj().T
 dist=float(np.sum(np.abs(np.linalg.eigvalsh(noisy-clean)))/2)
 assert dist<=3*gamma+1e-12
 outnoise.append({'gamma':gamma,'trace_distance':dist,'proved_cap':3*gamma})
result={'quantum_dynamics_oracles':records,'noise_oracles':outnoise,'source_dimension':len(labels),'H_nonzeros':H.nnz,
        'noiseless_numerical_propagation_method':'SciPy expm_multiply; ordinary numerical diagnostic, not rigorous acceptance',
        'seconds':time.monotonic()-start}
(ROOT/'development/numerical_validation.json').write_text(json.dumps(result,indent=2))
print('seconds',result['seconds'])
