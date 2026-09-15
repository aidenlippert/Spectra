"""Validation-only FCI upper witnesses for the straight-chain ladder."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,time,sys
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_implicit_certificate import rational_text

def interleaved(alpha,beta,n):
 s=sum(1<<(2*i) for i in range(n) if alpha&(1<<i))|sum(1<<(2*i+1) for i in range(n) if beta&(1<<i))
 c=sum(1 for i in range(n) if alpha&(1<<i) for j in range(n) if beta&(1<<j) and i>j)
 return s,(-1 if c%2 else 1)

def run(n, fixture_path, out, localized=False, top=1000):
 import numpy as np
 from pyscf import fci
 if top<1:raise ValueError('Positive top-k required')
 fixture=json.loads(Path(fixture_path).read_text()); t=time.monotonic()
 norb=fixture['modes']//2; ne=fixture['particles']; raw={tuple(tuple(x) for x in q['word']):F(q['coefficient']) for q in fixture['hamiltonian']}
 def coeff(word):
  from experiments.marginal_symbolic import canonical
  c=canonical({tuple(word):F(1)})
  return raw.get(next(iter(c)),F(0))/next(iter(c.values())) if c else F(0)
 h1=np.zeros((norb,norb)); eri=np.zeros((norb,)*4)
 for p in range(norb):
  for q in range(norb): h1[p,q]=float(coeff(((1,2*p),(0,2*q))))
  for r in range(norb):
   for q in range(norb):
    for s in range(norb): eri[p,r,q,s]=float(coeff(((1,2*p),(1,2*q+1),(0,2*s+1),(0,2*r))))
 e,ci=fci.direct_spin1.FCI().kernel(h1,eri,norb,(ne//2,ne//2),ecore=0.)
 strings=list(fci.cistring.make_strings(range(norb),ne//2)); vals=[]
 for i,a in enumerate(strings):
  for j,b in enumerate(strings):
   state,sign=interleaved(int(a),int(b),norb); vals.append((abs(float(ci[i,j])),state,int(round(float(ci[i,j])*sign*10**10))))
 vals.sort(reverse=True); full=len(vals); selected=vals if n<=6 else vals[:top]
 amps={s:a for _,s,a in selected if a}; states=sorted(amps); witness={'states':states,'amplitudes':[amps[s] for s in states]}
 if len(states)>4096:
  from research.certificate_scaling.streaming_reference_upper import upper as streaming_upper
  upper,upper_stats=streaming_upper(fixture,witness)
 else:
  oracle=DeterminantOracle(fixture); upper=oracle.upper(witness);upper_stats={'method':'existing_sparse_oracle'}
 hf=F(json.loads((Path(fixture_path).parent/'upper.json').read_text())['upper'])
 if len(selected)==full and abs(float(upper)-float(e))>1e-7: raise RuntimeError(f'FCI/CAR mapping mismatch: {float(upper)} vs {float(e)}')
 rec={'kind':'validation_only_fci_upper_v1','natoms':n,'localized':localized,'full_fci_support':full,'witness_support':len(states),'truncated':len(selected)<full,'selected_support':len(selected),'rounded_zero_count':len(selected)-len(states),'top_k':top if n>6 else None,'numeric_fci_electronic':float(e),'rational_witness_upper':rational_text(upper),'hf_upper':rational_text(hf),'hf_gap':float(hf-upper),'elapsed_seconds':time.monotonic()-t,'scope':'Validation upper only; prohibited as certificate-discovery input; top-k truncation and amplitude rounding are recorded separately.'}
 rec['upper_replay']=upper_stats
 out.mkdir(parents=True,exist_ok=True); (out/'upper.json').write_text(json.dumps({'independent_upper':witness,'upper':rational_text(upper)},indent=2)+'\n'); (out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n'); return rec

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--fixture-root',default='results/certificate_scaling/active_space_ladder'); ap.add_argument('--output-root',default='results/certificate_scaling/active_space_ladder_references_aligned'); ap.add_argument('--localized',action='store_true'); ap.add_argument('--atoms',type=int,nargs='+',default=[4,6,8,10]);ap.add_argument('--top',type=int,default=1000); a=ap.parse_args()
 root=Path(a.fixture_root); out=Path(a.output_root); allr=[]
 for n in a.atoms:
  allr.append(run(n,root/f'h{n}/fixture.json',out/f'h{n}',a.localized,a.top))
 (out/'summary.json').write_text(json.dumps(allr,indent=2)+'\n'); print(json.dumps(allr,indent=2))
if __name__=='__main__': main()
