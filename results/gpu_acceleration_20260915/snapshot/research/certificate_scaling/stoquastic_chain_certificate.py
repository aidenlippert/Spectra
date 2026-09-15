"""H-only positive-amplitude discovery and exact local-energy DP certificates.

Restricted spin chain backend, independent of the fermionic molecular solver.
No global state vector or configuration enumeration in discovery/verification.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,itertools,json,sys,time,hashlib

def index(a,b,c):return 4*a+2*b+c

def exponent(a,b,c):return 2*((a!=b)+(b!=c))-2

def parent_hamiltonian(n,hidden_q=F(2)):
    if n<4:raise ValueError('Periodic chain needs at least four sites')
    table=[str(hidden_q**exponent(a,b,c)) for a,b,c in itertools.product(range(2),repeat=3)]
    return {'sites':n,'diagonal_tables':[table[:] for _ in range(n)],'flip_rates':['1']*n}

def minmax_dp(tables):
    """Exact cyclic three-variable-factor DP. Returns extrema and state count."""
    n=len(tables);extrema=[];count=0
    for minimize in [True,False]:
        choose=min if minimize else max;final=[]
        for a,b in itertools.product(range(2),repeat=2):
            dp={(a,b):0}
            for j in range(2,n):
                new={}
                for (l,m),value in dp.items():
                    for r in range(2):
                        v=value+tables[j-1][index(l,m,r)];key=(m,r)
                        new[key]=v if key not in new else choose(new[key],v);count+=1
                dp=new
            final.extend(value+tables[n-1][index(l,m,a)]+tables[0][index(m,a,b)] for (l,m),value in dp.items())
        extrema.append(choose(final))
    return extrema[0],extrema[1],count

def local_energy(h,q,exact=True):
    cast=F if exact else lambda x: float(F(x))
    return [[cast(row[index(a,b,c)])-cast(t)*q**exponent(a,b,c) for a,b,c in itertools.product(range(2),repeat=3)]
             for row,t in zip(h['diagonal_tables'],h['flip_rates'])]

def verify(cert):
    start=time.monotonic();h=cert['hamiltonian'];n=h['sites'];q=F(cert['amplitude_domain_wall_weight'])
    if type(n) is not int or n<4 or q<=0 or len(h['diagonal_tables'])!=n or len(h['flip_rates'])!=n:raise ValueError('Invalid chain')
    if any(len(row)!=8 for row in h['diagonal_tables']):raise ValueError('Invalid local table')
    if any(F(t)<0 for t in h['flip_rates']):raise ValueError('Nonstoquastic off-diagonal rate')
    # Conditional 2x2 PSD flip blocks for every neighbor assignment.
    blocks=0
    for rate in h['flip_rates']:
        t=F(rate)
        for a,c in itertools.product(range(2),repeat=2):
            alpha=t*q**exponent(a,0,c);beta=t*q**exponent(a,1,c)
            assert alpha>=0 and beta>=0 and alpha*beta-t*t==0;blocks+=1
    tables=local_energy(h,q);lo,hi,count=minmax_dp(tables)
    # Weighted graph Laplacian plus diagonal local energy proves [lo,hi].
    return {'sites':n,'lower':str(lo),'upper':str(hi),'width':str(hi-lo),'lower_float':float(lo),'upper_float':float(hi),
            'width_float':float(hi-lo),'conditional_PSD_blocks_checked':blocks,'DP_transitions':count,
            'global_configurations_enumerated':0,'state_vector_constructed':False,'q':str(q),'replay_seconds':time.monotonic()-start}

def discover(h):
    from scipy.optimize import minimize_scalar
    import math
    start=time.monotonic();calls=0;transitions=0
    def objective(theta):
        nonlocal calls,transitions
        lo,_,count=minmax_dp(local_energy(h,math.exp(theta),False));calls+=1;transitions+=count
        return -lo
    result=minimize_scalar(objective,bounds=(-math.log(4),math.log(4)),method='bounded',options={'xatol':1e-12,'maxiter':200})
    q=F(round(math.exp(result.x)*10**6),10**6)
    cert={'kind':'stoquastic_periodic_positive_factor_v1','hamiltonian':h,'amplitude_domain_wall_weight':str(q)}
    rec=verify(cert);rec.update({'discovery_seconds':time.monotonic()-start,'objective_calls':calls,'discovery_DP_transitions':transitions,
             'optimizer_success':bool(result.success),'parameter_box_q':['1/4','4'],'source_state_or_weights_used':False,
             'scope':'Uniform positive Jastrow family on stoquastic periodic spin chains; not generic fermions.'})
    return cert,rec

def control():
    h=parent_hamiltonian(5)
    for q in [F(2),F(3,2),F(5,3)]:
        tables=local_energy(h,q);lo,hi,_=minmax_dp(tables)
        brute=[sum(tables[i][index(x[i-1],x[i],x[(i+1)%5])] for i in range(5)) for x in itertools.product(range(2),repeat=5)]
        assert (lo,hi)==(min(brute),max(brute))
    cert={'hamiltonian':h,'amplitude_domain_wall_weight':'2'};assert verify(cert)['width']=='0'
    bad=dict(h,flip_rates=['-1']+h['flip_rates'][1:])
    try:verify(dict(cert,hamiltonian=bad))
    except ValueError:pass
    else:raise AssertionError('Nonstoquastic input accepted')
    return {'DP_bruteforce_small_control':True,'nonstoquastic_rejected':True}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',type=Path);p.add_argument('--out',type=Path,default=Path('results/certificate_scaling/stoquastic_chain'));a=p.parse_args()
    if a.verify:print(json.dumps(verify(json.loads(a.verify.read_text())),indent=2))
    else:
        controls=control();a.out.mkdir(parents=True,exist_ok=True);rows=[]
        for n in [8,16,32,64,128,256]:
            h=parent_hamiltonian(n);cert,rec=discover(h);f=a.out/f'L{n}.json';f.write_text(json.dumps(cert,separators=(',',':'))+'\n')
            rec.update({'certificate_bytes':f.stat().st_size,'certificate_sha256':hashlib.sha256(f.read_bytes()).hexdigest()});rows.append(rec)
            print(json.dumps(rec),flush=True)
        (a.out/'summary.json').write_text(json.dumps({'controls':controls,'results':rows},indent=2)+'\n')
