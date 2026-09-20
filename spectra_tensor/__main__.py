"""Command-line entry point; verification works under python -S."""
from pathlib import Path
import argparse,json,sys

def main():
    if len(sys.argv)>1 and sys.argv[1]=='verify':
        sys.argv.pop(1)
        from .exact import main as verify_main
        verify_main();return
    p=argparse.ArgumentParser(description='Construct and verify finite Hubbard responses without sector enumeration.')
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('request');r.add_argument('--sites',type=int,required=True);r.add_argument('--U',default='8');r.add_argument('--leg',default='1');r.add_argument('--rung',default='1');r.add_argument('--omega',nargs='+',default=['2','-1','8']);r.add_argument('--eta',default='3/2');r.add_argument('--target',default='1/1000');r.add_argument('--out',type=Path,required=True)
    s=sub.add_parser('solve');s.add_argument('--request',type=Path,required=True);s.add_argument('--out',type=Path,required=True);s.add_argument('--bonds',nargs='+',type=int,default=[48,96,160]);s.add_argument('--seed-steps',type=int,default=12);s.add_argument('--sweeps',type=int,default=4);s.add_argument('--backend',choices=['gmp','python'],default='gmp')
    sub.add_parser('verify',help='Independent replay; use verify --help')
    a=p.parse_args()
    if a.command=='request':
        from .exact import rat,validate_model,source
        n=a.sites
        if n%2 or not 2<=n<=200:p.error('even sites in [2,200] required')
        spec=dict(kind='spin_independent_hubbard_v1',sites=n,particles=n,U=[str(rat(a.U))]*n,edges=[[2*j,2*j+1,str(rat(a.rung))] for j in range(n//2)]+[[2*j+s,2*j+2+s,str(rat(a.leg))] for j in range(n//2-1) for s in (0,1)],mode_order='up_then_down',units='t',boundary='open')
        validate_model(spec)
        if rat(a.eta)<=0 or rat(a.target)<=0:p.error('positive eta and target required')
        request=dict(kind='tensor_hubbard_response_request_v1',model=spec,source=source(spec),frequencies=[dict(omega=str(rat(w)),eta=str(rat(a.eta))) for w in a.omega],target_radius=str(rat(a.target)))
        a.out.parent.mkdir(parents=True,exist_ok=True)
        with a.out.open('x') as f:json.dump(request,f,indent=2)
        print(a.out)
    elif a.command=='solve':
        from .api import solve
        result=solve(json.loads(a.request.read_text()),a.out,bonds=tuple(a.bonds),seed_steps=a.seed_steps,sweeps=a.sweeps,backend=a.backend)
        print(json.dumps({k:v for k,v in result.items() if k!='queries'},indent=2))
        if result['status']!='target_met':sys.exit(2)
if __name__=='__main__':main()
