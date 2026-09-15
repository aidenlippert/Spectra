"""Simple HF-occupation parent test; no residual-engineered parent."""
import importlib.util, json, numpy as np, math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path('/Users/aidenlippert/Documents/Spectra')
spec=importlib.util.spec_from_file_location('o',ROOT/'research/all_angles_20260913/frontier/exact_h4_oracle.py')
o=importlib.util.module_from_spec(spec); spec.loader.exec_module(o)
def main():
 d=json.loads((ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json').read_text()); A,states=o.matrix(d); n=len(A)
 hf=min(range(n),key=lambda i:A[i][i]); E=A[hf][hf]
 P=[[F(int(i==j and i!=hf)) for j in range(n)] for i in range(n)]
 R=[[A[i][j]-(E if i==j else 0)-P[i][j] for j in range(n)] for i in range(n)]
 vals=np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in A])); rows=[]
 for eta in (F(0),F(1,2),F(1)):
  C=np.array([[float(R[i][j]+eta*P[i][j]) for j in range(n)] for i in range(n)])
  eps=F(max(0,math.ceil(-np.linalg.eigvalsh(C)[0]*10**6)),10**6)
  while eps and o.ldl_psd([[R[i][j]+eta*P[i][j]+(eps if i==j else 0) for j in range(n)] for i in range(n)])[0]: eps-=F(1,10**6)
  while not o.ldl_psd([[R[i][j]+eta*P[i][j]+(eps if i==j else 0) for j in range(n)] for i in range(n)])[0]: eps+=F(1,10**6)
  rows.append({'eta':str(eta),'epsilon':str(eps),'width':str(eps),'lower':str(E-eps)})
 out={'dimension':n,'hf_state':states[hf],'parent':'global projector I-|HF><HF| (diagnostic, not local sum)','E_ref':str(E),'results':rows,'true_ground_numeric':float(vals[0]),'interpretation':'exact relative certificates; widths are eps'}
 (ROOT/'results/wave2_20260913/theory/parent_test.json').write_text(json.dumps(out,indent=2)); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
