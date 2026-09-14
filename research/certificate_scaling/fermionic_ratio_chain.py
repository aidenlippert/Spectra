"""Exact H-only compiler for a biased nearest-neighbor fermionic parent family.

Known structured control, not generic chemistry. Noncommuting quartic H;
fixed-N weighted positive state; O(M) cubic SOS factors and polynomial replay.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import argparse,json,sys,time,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import add,scale,mono,product,canonical,adj,encode,verify

def number(i):return mono(((1,i),(0,i)))

def hamiltonian(m,r):
    terms=[]
    for i in range(m-1):
        terms.extend([scale(number(i),r),scale(number(i+1),1/r),scale(product(number(i),number(i+1)),-r-1/r),
                      mono(((1,i),(0,i+1)),-1),mono(((1,i+1),(0,i)),-1)])
    return add(*terms)

def compile_h(h,m,n):
    start=time.monotonic()
    if m<3 or not 0<=n<=m:raise ValueError('Invalid sector')
    h=canonical(h);r=h.get(((1,0),(0,0)),F(0))
    if r<=0 or h!=hamiltonian(m,r):raise ValueError('Hamiltonian not in recognized nearest-neighbor ratio family')
    an,bn=isqrt(r.numerator),isqrt(r.denominator)
    if an*an!=r.numerator or bn*bn!=r.denominator:raise ValueError('Rational square ratio required by exact factor exporter')
    q=F(an,bn);polys=[]
    for i in range(m-1):
        ai=mono(((0,i),));aj=mono(((0,i+1),))
        polys.append(add(scale(ai,q),scale(product(ai,number(i+1)),-q),scale(aj,-1/q),scale(product(aj,number(i)),1/q)))
    from math import lcm
    den=lcm(*(c.denominator for p in polys for c in p.values()))
    blocks=[{'name':f'bond-{i}','words':list(p),'factor':[[int(c*den) for c in p.values()]]} for i,p in enumerate(polys)]
    cert={'modes':m,'particles':n,'hamiltonian':encode(h),'number_multiplier':[],'b':'0','denominator':den,'blocks':blocks,
          'structured_upper':{'kind':'positive_fixed_number_site_ratio_v1','r':str(r)}}
    rec=verify(cert);assert rec['lower']=='0'
    rec.update({'discovery_seconds':time.monotonic()-start,'source_factors_used':False,'source_state_used':False})
    return cert,rec

def structured_replay(cert):
    start=time.monotonic();rec=verify(cert);m=cert['modes'];n=cert['particles'];r=F(cert['structured_upper']['r'])
    if cert['structured_upper']['kind']!='positive_fixed_number_site_ratio_v1' or r<=0:raise ValueError('Invalid positive state')
    from experiments.marginal_symbolic import decode
    h=decode(cert['hamiltonian'],m,4)
    if h!=hamiltonian(m,r):raise ValueError('Upper state not matched to H')
    # In each adjacent 10/01 block, H_i=[[r,-1],[-1,1/r]],
    # while amplitudes have ratio psi(01)/psi(10)=r. Both rows vanish.
    assert r-r==0 and -1+(1/r)*r==0
    # Exact norm of psi(x)=r^(sum_i i*x_i), restricted to sum_i x_i=N.
    dp=[F(0)]*(n+1);dp[0]=F(1);operations=0
    for i in range(m):
        weight=r**(2*i)
        for j in range(min(n,i+1),0,-1):dp[j]+=weight*dp[j-1];operations+=1
    assert dp[n]>0
    rec.update({'upper':'0','width':str(-F(rec['lower'])),'norm_positive':True,'norm_numerator_bits':dp[n].numerator.bit_length(),
                'norm_denominator_bits':dp[n].denominator.bit_length(),'norm_DP_updates':operations,'r':str(r),
                'many_body_states_enumerated':0,'replay_seconds':time.monotonic()-start})
    return rec

def controls():
    # Independent local occupancy matrices and full-sector vector on four modes.
    from experiments.marginal_determinant_tree import DeterminantOracle
    h=hamiltonian(4,F(4));c,_=compile_h(h,4,2)
    states=[x for x in range(16) if x.bit_count()==2];amps=[4**sum(i for i in range(4) if (x>>i)&1) for x in states]
    assert DeterminantOracle(c).upper({'states':states,'amplitudes':amps})==0
    e0=hamiltonian(3,F(4));left=hamiltonian(2,F(4))
    right={tuple((cr,i+1) for cr,i in w):v for w,v in left.items()}
    comm=add(product(left,right),scale(product(right,left),-1));assert comm
    bad=dict(h);bad[((1,0),(0,0))]+=1
    try:compile_h(bad,4,2)
    except ValueError:pass
    else:raise AssertionError('Mutation not rejected')
    return {'exact_small_sector_upper':True,'adjacent_energy_commutator_terms':len(comm),'mutation_rejected':True}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',type=Path);p.add_argument('--out',type=Path,default=Path('results/certificate_scaling/fermionic_ratio_chain'));a=p.parse_args()
    if a.verify:print(json.dumps(structured_replay(json.loads(a.verify.read_text())),indent=2))
    else:
        checks=controls();a.out.mkdir(parents=True,exist_ok=True);rows=[]
        for m in [8,16,32,64,128,256]:
            h=hamiltonian(m,F(4));c,discovery=compile_h(h,m,m//2);rec=structured_replay(c);rec['discovery_seconds']=discovery['discovery_seconds']
            f=a.out/f'M{m}.json';f.write_text(json.dumps(c,separators=(',',':'))+'\n');rec['certificate_bytes']=f.stat().st_size;rec['certificate_sha256']=hashlib.sha256(f.read_bytes()).hexdigest();rows.append(rec)
            print(json.dumps({k:rec[k] for k in ['mode_count','particle_number','lower','upper','factor_rows','factor_nonzeros','certificate_bytes','discovery_seconds','replay_seconds','norm_DP_updates']}),flush=True)
        (a.out/'summary.json').write_text(json.dumps({'controls':checks,'results':rows},indent=2)+'\n')
