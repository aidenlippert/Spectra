from pathlib import Path
import argparse,json,resource,sys
from .reference import solve_reference

def main():
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    rows=[]
    for n in [8,12,20,50]:
        spec=dict(kind='spin_independent_hubbard_v1',sites=n,particles=n,U=['8']*n,edges=[[2*j,2*j+1,'1'] for j in range(n//2)]+[[2*j+s,2*j+2+s,'1'] for j in range(n//2-1) for s in (0,1)],mode_order='up_then_down',units='t',boundary='open')
        row=dict(sites=n,**solve_reference(spec,2+1.5j));rows.append(row);print(row,flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x') as f:json.dump(dict(results=rows,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),timing_note='Development runs may overlap; no isolated speedup claim.'),f,indent=2)
if __name__=='__main__':main()
