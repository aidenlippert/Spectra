"""Independent ordinary-CAR comparison of every H8 full moment through14."""
from pathlib import Path
import json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
ROOT=Path(__file__).resolve().parents[1]
def main():
    started=time.monotonic();h=json.loads((ROOT/'hamiltonian.json').read_text());c=json.loads((ROOT/'singlet_valence_embedding/certificate.json').read_text())
    V=[{int(s):a for s,a in v.items()} for v in c['basis']];oracle=DeterminantOracle(h)
    def action(v):
        out={}
        for s,a in v.items():
            for t,b in oracle.action(s).items():
                if b.denominator!=1:raise ValueError('Integer model expected')
                out[t]=out.get(t,0)+a*b.numerator
        return {s:a for s,a in out.items() if a}
    powers=[V]
    for _ in range(7):powers.append([action(v) for v in powers[-1]])
    M=[]
    for k in range(15):
        left,right=powers[k//2],powers[k-k//2]
        M.append([[sum(a*v.get(s,0) for s,a in u.items()) for v in right] for u in left])
    quotient=SymmetryMomentOracle(h);N,r=quotient.moments(V,14)
    if M!=N:raise ValueError('H8 full moments differ from original CAR')
    receipt={'accepted':True,'orders':[0,14],'all_2940_entries_equal':True,
             'original_car_source_actions':len(oracle.cache),'original_car_referenced':oracle.referenced_state_count(),
             'quotient_work':r,'seconds':time.monotonic()-started,
             'scope':'Independent finite H8 ordinary-CAR equality check for every full moment through14; production moment replay itself uses the signed-orbit route.'}
    (ROOT/'symmetry_moments/independent_car_comparison.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
