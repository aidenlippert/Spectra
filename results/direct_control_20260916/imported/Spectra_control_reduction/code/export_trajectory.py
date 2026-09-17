from propose import *
import hashlib

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def export_case(name,basis_file,d,phases,robust_amp='0',trace_radius='0',threshold='-3/5'):
 t=time.monotonic();labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');P=np.load(ROOT/'development'/basis_file)['J'][:,:d]
 JD=2**32; HD=2**40;ZD=2**48
 Jint=np.rint(P*JD).astype(np.int64);J=Jint.astype(float)/JD
 hs=[J.T@K@J,J.T@(D[:,None]*J),J.T@W@J]
 # Exact model matrices are proposal aids only; literal residual acceptance proves validity.
 z=J.T@psi; degree=18; dt=.25; segments=[];hcache={}
 for duration,u,v in phases:
  key=(str(F(u)),str(F(v)))
  if key not in hcache:
   hn=np.rint((hs[0]+float(F(u))*hs[1]+float(F(v))*hs[2])*HD).astype(np.int64)
   hn=np.triu(hn)+np.triu(hn,1).T
   hcache[key]=hn
  hi=hcache[key];hh=hi.astype(float)/HD;ev,V=linalg.eigh(hh)
  steps=round(duration/dt);assert steps*dt==duration
  for step in range(steps):
   coeff=[z.copy()]
   for l in range(1,degree+1):coeff.append((-1j*dt/l)*(hh@coeff[-1]))
   cre=np.rint(np.array(coeff).real*ZD).astype(np.int64).tolist();cim=np.rint(np.array(coeff).imag*ZD).astype(np.int64).tolist()
   segments.append({'duration':'1/4','u':key[0],'v':key[1],'re':cre,'im':cim})
   z=V@(np.exp(-1j*dt*ev)*(V.T@z))
 data={'kind':'fixed_number_control_trajectory_v1','modes':16,'spin_counts':[4,4],
  'fixture_sha256':sha(ROOT/'inputs/fixture.json'),'state_sha256':sha(ROOT/'inputs/state.json'),
  'basis_denominator':JD,'basis':Jint.tolist(),'model_denominator':HD,'models':[{'u':k[0],'v':k[1],'h':a.tolist()} for k,a in hcache.items()],
  'trajectory_denominator':ZD,'segments':segments,'energy_shift':'1851/200',
  'control_operator_norm_bounds':['2','2'],'robust_amplitude_Ha':robust_amp,'initial_trace_distance':trace_radius,
  'target_D_upper':threshold,'target_D_lower':'1' if name.startswith('baseline') else '-2',
  'integrated_dephasing_rate':'0' if name.startswith('baseline') else '1/1000',
  'scope':'Actual rational H8 driven from its specified rational MPS; finite-model controls, not a laboratory pulse or true-ground initialization.',
  'discovery':'Enumerated 4900 balanced-spin determinants and built numerical training snapshots; no enumeration-free discovery or scaling claim.'}
 out=ROOT/'inputs'/f'{name}.json';out.write_text(json.dumps(data,separators=(',',':'))+'\n')
 print(name,'d',d,'segments',len(segments),'bytes',out.stat().st_size,'seconds',time.monotonic()-t,flush=True)
if __name__=='__main__':
 cs=np.load(ROOT/'development/adapted_basis.npz')['controls'];phases=[(.5,F(int(round(u*1e6)),10**6),F(int(round(v*1e6)),10**6)) for u,v in cs]
 export_case('short24','adapted_basis.npz',24,phases,'1/1000','1/2000')
 export_case('short8_refusal','adapted_basis.npz',8,phases,'1/1000','1/2000')
 export_case('long64','long_adapted_basis.npz',64,[(21,F(1,5),F(1,10))],'9/100000','1/2000')
 export_case('baseline32','baseline_basis.npz',32,[(21,F(0),F(0))],'0','0','2')
