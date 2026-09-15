from fractions import Fraction as F
import json, numpy as np
from research.all_angles_20260913.frontier.exact_h4_oracle import matrix, ldl_psd

def main():
 data=json.load(open('results/certificate_scaling/active_space_ladder/h4/fixture.json')); H,_=matrix(data); n=70
 D=[[H[i][j] if i==j else F(0) for j in range(n)] for i in range(n)]; V=[[H[i][j]-D[i][j] for j in range(n)] for i in range(n)]
 rows=[]
 for lam in [F(0),F(1,4),F(1,2),F(3,4),F(1)]:
  lo=min(D[i][i]+lam*V[i][i]-lam*sum(abs(V[i][j]) for j in range(n) if j!=i) for i in range(n)); up=D[0][0]+lam*V[0][0]
  A=np.array([[float(D[i][j]+lam*V[i][j]) for j in range(n)] for i in range(n)])
  rows.append({'lambda':str(lam),'gershgorin_lower':str(lo),'trial_e0_upper':str(up),'certified_width':float(up-lo),'oracle_ground':float(np.linalg.eigvalsh(A)[0])})
 dmin=min(D[i][i] for i in range(n)); P=[[D[i][j]-(dmin if i==j else 0) for j in range(n)] for i in range(n)]; rel=[]
 for eta in [F(0),F(1,10),F(1,2),F(1)]:
  for eps in [F(1,100),F(1,10),F(1)]:
   S=[[V[i][j]+eta*P[i][j]+eps*(i==j) for j in range(n)] for i in range(n)]; ok,piv,idx=ldl_psd(S); rel.append({'eta':str(eta),'eps':str(eps),'status':'PASS' if ok else 'FAIL','pivot':str(idx if idx is not None else piv)})
 out={'fixture_sha256':'8e64505deb4cc06a87f75e669b0f3073d15ebca72daa6d086e06eeea0b0f5120','dimension':n,'path':rows,'relative_scan':rel,'V_max_abs_row_sum':str(max(sum(abs(x) for x in r) for r in V)),'zero_pivot_regression':ldl_psd([[F(0),F(1)],[F(1),F(0)]])[0]}; json.dump(out,open('results/wave2_20260913/continuation/receipt.json','w'),indent=2); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
