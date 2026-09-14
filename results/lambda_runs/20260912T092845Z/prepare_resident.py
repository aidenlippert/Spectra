from pathlib import Path
import json,time,sys
from fractions import Fraction as F
import numpy as np
from results.marginal_graded_hubbard8.discovery import joint_gpu_benchmark as b
out=Path(sys.argv[1]);start=time.perf_counter()
cert=json.loads((b.BASE/'six_site_projector/refined_certificate.json').read_text());local=cert['local_window']
origin=[F(local[k][j]) for k,j in [('onsite_profile',0),('onsite_profile',1),('hopping_profile',0),('hopping_profile',1),('density_profile',0),('density_profile',1)]]
half,hn=b._projector_vector(cert['vector'],6)
cs=json.loads((b.BASE/'charged_projector/source.json').read_text())['physical_state']
cf,cn=b.charged_vectors(cs)
mats,der=b._profiles_resident(origin,local,half,hn,cf,cn,F(1,1000))
grams=b._gram_batch(cert['vector'],cs,0.5)
arrays={'origin':np.array(origin,dtype=float),'alpha':np.array(float(F(cert['penalty']))),'beta':np.array(.001),'ratio':np.array(.5)}
for i,(_,a,p,q) in enumerate(mats):
 arrays[f'A{i}']=a;arrays[f'P{i}']=p;arrays[f'Q{i}']=q;arrays[f'D{i}']=np.array([der[j][i] for j in range(6)])
for i,g in enumerate(grams):arrays[f'G{i}']=g
# Independently reconstruct changed physical matrices to check the derivative assembly.
max_error=0.
for offsets in [[F(j-2,1000) for j in range(6)],[F(2-j,2000) for j in range(6)]]:
 changed=b._resident([v+d for v,d in zip(origin,offsets)],local,half,hn,cf,cn)
 for i,(_,a,p,q) in enumerate(changed):
  assembled=mats[i][1]+np.einsum('i,ijk->jk',np.array(offsets,dtype=float),arrays[f'D{i}'])
  max_error=max(max_error,float(np.max(np.abs(a-assembled))))
assert max_error<1e-10,max_error
np.savez(out/'resident.npz',**arrays)
receipt={'preparation_seconds':time.perf_counter()-start,'sectors':len(mats),'gram_blocks':len(grams),'matrix_bytes':sum(a.nbytes for a in arrays.values()),'direct_assembly_max_abs_error':max_error,'shapes':{str(n):sum(len(a)==n for _,a,_,_ in mats) for n in sorted({len(a) for _,a,_,_ in mats})}}
(out/'resident_preparation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
