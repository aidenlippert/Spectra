"""Exact small-chain validation of an ordinary occupation-basis norm transfer."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from tests.test_marginal_hopping_filter import action,hopping,dot
from experiments.marginal_local_hubbard_block import _actions

def main():
    # This diagnostic intentionally reuses the independent direct-CAR test
    # helpers. It is not a production verifier or an H8 energy certificate.
    phi={3:1,12:1,9:2,6:-2}
    h=_actions(2,0,1)
    def move(v):
        out={}
        for s,a in v.items():
            for t,b in h[s].items():out[t]=out.get(t,0)+a*b
        return out
    # Pair index is ket+4*bra. A two-site block has no middle to trace.
    G=[[phi.get(l+4*r,0)*phi.get(lp+4*rp,0) for rp in range(4) for r in range(4)]
       for lp in range(4) for l in range(4)]
    e=[int(i%4==i//4) for i in range(16)]
    def vm(v,M):return [sum(v[i]*M[i][j] for i in range(16)) for j in range(16)]
    cases=[]
    for a,b in [(F(1,5),0),(F(1,5),F(1,7)),(0,1)]:
        gates={}
        for s in range(16):
            h1=h[s];h2=move(h1)
            gates[s]={t:int(t==s)-a*h1.get(t,0)+b*h2.get(t,0) for t in range(16)}
        f2={s:{t:sum(gates[s][k]*gates[t][k] for k in range(16)) for t in range(16)} for s in range(16)}
        B=[[f2[r+4*l][rp+4*lp] for lp in range(4) for l in range(4)] for rp in range(4) for r in range(4)]
        row=vm(e,G)
        for q in range(1,5):
            v={0:1}
            for j in range(q):
                v={s+(t<<(4*j)):x*y for s,x in v.items() for t,y in phi.items()}
            for cut in range(2,2*q,2):
                h1=action(hopping(cut),v);h2=action(hopping(cut),h1)
                v={s:v.get(s,0)-a*h1.get(s,0)+b*h2.get(s,0) for s in v.keys()|h1.keys()|h2.keys()}
            direct=dot(v,v);transfer=sum(x*y for x,y in zip(row,e))
            if direct!=transfer:raise ValueError('Transfer norm disagrees with direct CAR')
            cases.append({'a':str(a),'b':str(b),'blocks':q,'norm':str(transfer)})
            row=vm(vm(row,B),G)
    sources=['results/marginal_graded_hubbard8/discovery/boundary_norm_transfer_probe.py',
        'tests/test_marginal_hopping_filter.py','experiments/marginal_local_hubbard_block.py',
        'experiments/marginal_transfer_verify.py']
    result={'accepted':True,'dimension':16,'cases':cases,
        'formula':'e^T G (B G)^(q-1) e; e_(l,lprime)=delta_(l,lprime)',
        'scope':'Exact norm-only agreement for a specified two-site block, three gates, and one through four blocks. No energy-insertion, general accuracy/cost, or large-chain interval certificate yet.',
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    path=ROOT/'results/marginal_graded_hubbard8/adapted_block/boundary_norm_transfer_probe.json'
    path.write_text(json.dumps(result,indent=2)+'\n');print('Accepted12 exact norm comparisons',flush=True)

if __name__=='__main__':main()
