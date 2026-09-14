"""Fixed sparse-action-budget probes; these are single-boundary moments only."""
from itertools import product
from pathlib import Path
import json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
ROOT=Path(__file__).resolve().parents[1]
def boundary(sites):
    v={}
    for choices in product((0,1),repeat=sites//2):
        s=0
        for pair,c in enumerate(choices):s|=1<<(4*pair+c);s|=1<<(4*pair+2+1-c)
        v[s]=(-1)**sum(choices)
    return v
def main():
    out=ROOT/'symmetry_moment_scaling';out.mkdir(exist_ok=True);rows=[]
    for sites in [8,12,16]:
        o=SymmetryMomentOracle(build(sites,4,1));v=boundary(sites)
        for order in [4,8,12,18,24]:
            start=time.monotonic()
            try:M,r=o.moments([v],order)
            except ValueError as error:
                rows.append({'sites':sites,'moment_order':order,'accepted':False,'reason':str(error),'source_actions':len(o.base.cache),'seconds':time.monotonic()-start});break
            else:rows.append({'sites':sites,'moment_order':order,'accepted':True,'boundary_determinants':len(v),'last_moment':str(M[-1][0][0]),'work':r,'seconds':time.monotonic()-start})
            finally:(out/'receipt.json').write_text(json.dumps({'results':rows,'scope':'Single product-of-singlets boundary moment probes at varying order under the unchanged4096 determinant-action cap. No larger-system energy certificate or equal-accuracy scaling claim.'},indent=2)+'\n')
            print(rows[-1],flush=True)
if __name__=='__main__':main()
