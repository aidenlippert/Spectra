"""Standard-library replay of the one-block integer congruence witness."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_charged_projectors import joint_overlap_grams
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'joint_projector/congruence_probe'


def main():
    started=time.monotonic();proposal_path=OUT/'one_block_receipt.json';proposal=json.loads(proposal_path.read_text())
    source_path=BASE/'joint_projector/signed_density/profile_joint_r1_2_certificate.json';c=json.loads(source_path.read_text())
    data=joint_overlap_grams(c['vector'],c['joint']['vector'],4)[tuple(proposal['sector'])]
    r=F(c['joint']['ratio']);B=F(c['joint']['projector_sum_ceiling']);n=len(data['gram'])
    A=[[(B*data['norms'][i]/(1 if data['sources'][i]==0 else r) if i==j else 0)-v for j,v in enumerate(row)] for i,row in enumerate(data['gram'])]
    scale=lcm(*(F(v).denominator for row in A for v in row));A=[[int(v*scale) for v in row] for row in A]
    if hashlib.sha256(json.dumps(A,separators=(',',':')).encode()).hexdigest()!=proposal['matrix_sha256']:raise ValueError('Fresh Gram matrix differs from proposed matrix')
    R=proposal['triangular_integer_congruence']
    if len(R)!=n or any(len(row)!=n for row in R) or any(type(v) is not int for row in R for v in row):raise ValueError('Integer square congruence required')
    if any(R[i][j] for i in range(n) for j in range(i)) or any(not R[i][i] for i in range(n)):raise ValueError('Nonsingular triangular congruence required')
    # Generic dot products, independently of the triangular loop optimization
    # used during discovery. Positive common denominators do not affect PSD.
    Rt=list(zip(*R));AR=[[sum(a*b for a,b in zip(row,col)) for col in Rt] for row in A]
    ARt=list(zip(*AR));C=[[sum(a*b for a,b in zip(row,col)) for col in ARt] for row in Rt]
    if any(C[i][j]!=C[j][i] for i in range(n) for j in range(n)):raise ValueError('Congruence symmetry failed')
    margins=[row[i]-sum(abs(v) for j,v in enumerate(row) if i!=j) for i,row in enumerate(C)]
    if min(margins)<=0:raise ValueError('Strict diagonal dominance failed')
    receipt={'accepted':True,'dimension':n,'sector':proposal['sector'],'smallest_exact_margin':str(min(margins)),
             'seconds':time.monotonic()-started,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),proposal_path,source_path,ROOT/'experiments/marginal_charged_projectors.py')},
             'scope':'Independent standard-library reconstruction and generic integer congruence multiply verify one actual joint Gram block. Not an entire energy replay or a singular-PSD verifier.'}
    (OUT/'one_block_stdlib_replay.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
