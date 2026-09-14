"""Bounded exact diagonal-dominance congruence probe on one actual Gram block.

NumPy/SciPy propose a rational triangular change of basis. Only exact integer
arithmetic establishes strict diagonal dominance of the congruent matrix.
"""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import hashlib,json,sys,time
import numpy as np
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_charged_projectors import joint_overlap_grams
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'joint_projector/congruence_probe'


def main():
    OUT.mkdir(exist_ok=True);started=time.monotonic()
    c=json.loads((BASE/'joint_projector/signed_density/profile_joint_r1_2_certificate.json').read_text())
    grams=joint_overlap_grams(c['vector'],c['joint']['vector'],4);key,data=max(grams.items(),key=lambda item:len(item[1]['gram']))
    r=F(c['joint']['ratio']);B=F(c['joint']['projector_sum_ceiling']);n=len(data['gram'])
    A=[[(B*data['norms'][i]/(1 if data['sources'][i]==0 else r) if i==j else 0)-v for j,v in enumerate(row)] for i,row in enumerate(data['gram'])]
    denominator=lcm(*(F(v).denominator for row in A for v in row));ints=[[int(v*denominator) for v in row] for row in A]
    preparation=time.monotonic()-started;prop=time.monotonic();scale=max(ints[i][i] for i in range(n))
    numeric=np.array(ints,float)/float(scale);chol=np.linalg.cholesky(numeric)
    R=solve_triangular(chol.T,np.eye(n),lower=False);Rint=[[round(float(v)*10**12) for v in row] for row in R]
    if any(Rint[i][j] for i in range(n) for j in range(i)) or any(Rint[i][i]==0 for i in range(n)):raise ValueError('Singular or nontriangular proposed congruence')
    proposal=time.monotonic()-prop;exact_start=time.monotonic()
    # B=R^T A R; common positive denominator may be omitted for PSD.
    AR=[[sum(ints[i][k]*Rint[k][j] for k in range(j+1)) for j in range(n)] for i in range(n)]
    congruent=[[sum(Rint[k][i]*AR[k][j] for k in range(i+1)) for j in range(n)] for i in range(n)]
    if any(congruent[i][j]!=congruent[j][i] for i in range(n) for j in range(n)):raise ValueError('Exact congruence lost symmetry')
    margins=[row[i]-sum(abs(v) for j,v in enumerate(row) if j!=i) for i,row in enumerate(congruent)]
    if min(margins)<=0:raise ValueError('Exact diagonal dominance refused proposed congruence')
    receipt={'accepted':True,'sector':key,'dimension':n,'matrix_scale':str(denominator),'matrix_sha256':hashlib.sha256(json.dumps(ints,separators=(',',':')).encode()).hexdigest(),'triangular_integer_congruence':Rint,'smallest_exact_dominance_margin':str(min(margins)),'preparation_seconds':preparation,'numeric_proposal_seconds':proposal,'exact_integer_seconds':time.monotonic()-exact_start,'scope':'One actual tight r=1/2 Gram sector only. Nonsingular upper-triangular integer congruence and exact strict diagonal dominance prove this block positive definite. No other sector, whole energy certificate, singular PSD case or controlled speedup benchmark is covered.'}
    (OUT/'one_block_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='triangular_integer_congruence'}),flush=True)


if __name__=='__main__':main()
