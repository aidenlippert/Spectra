"""Exact fixed-N Gershgorin residual bounds on exterior-power coefficients.

This strengthens certificate replay; it does not infer small certificate
dictionaries or replace the cost of finding the original SOS factors.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import verified_residual


def residual_bounds(residual,m,n):
    bodies={}
    for word,c in residual.items():
        if not word:continue
        k=len(word)//2
        if len(word)!=2*k or [cr for cr,_ in word]!=[1]*k+[0]*k:
            raise ValueError('Balanced normal-ordered residual required')
        i=tuple(p for _,p in word[:k]);j=tuple(p for _,p in word[k:])
        bodies.setdefault(k,{})[i,j]=c*((-1)**(k*(k-1)//2))
    total=F(residual.get((),0));details=[]
    for k,matrix in sorted(bodies.items()):
        if any(matrix.get((j,i),F(0))!=c for (i,j),c in matrix.items()):
            raise ValueError('Non-Hermitian residual coefficient block')
        norms={};diagonal={}
        for (i,j),c in matrix.items():
            norms.setdefault(i,F(0))
            if i==j:diagonal[i]=c
            else:norms[i]+=abs(c)
        endpoints=[diagonal.get(i,F(0))-row_sum for i,row_sum in norms.items()]
        dimension=comb(m,k)
        if len(norms)<dimension:endpoints.append(F(0))
        lower=min(endpoints,default=F(0));l1=sum(abs(c) for c in matrix.values())
        sector_lower=lower*(comb(n,k) if k<=n else 0)
        chosen=max(-l1,sector_lower);total+=chosen
        details.append({'body':k,'wedge_dimension':dimension,'active_rows':len(norms),
                        'coefficient_l1':str(l1),'gershgorin_lower':str(lower),
                        'fixed_N_lower':str(sector_lower),'chosen_lower':str(chosen)})
    return total,details


def replay(cert):
    start=time.monotonic();residual,old=verified_residual(cert);m,n=cert['modes'],cert['particles']
    if 'number_multiplier' not in cert:raise ValueError('This replay requires explicit multiplier coefficients')
    correction,details=residual_bounds(residual,m,n)
    lower=F(cert['b'])+correction
    if lower<F(old['lower']):raise AssertionError('Alternative bound weaker than coefficient L1')
    return {'method':'fixed_number_wedge_gershgorin_v1','modes':m,'particles':n,
            'coefficient_lower':old['lower'],'lower':str(lower),'lower_float':float(lower),
            'improvement':str(lower-F(old['lower'])),'constant_residual':str(residual.get((),0)),
            'bodies':details,'replay_seconds':time.monotonic()-start,'many_body_states_enumerated':0,
            'scope':'fixed total N, no spin-sector assumption; original SOS plus exact residual inequalities'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    rec=replay(json.loads(a.certificate.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
