"""Read-only rational validation of reduced driven-H8 trajectories.

No numerical libraries. The validation deliberately enumerates the 4900-state
balanced-spin sector and charges that dependency. It does not claim scalable
many-body discovery. The final finite-time and waveform-neighborhood bounds
follow from unitary defect propagation, not sampled dynamical errors.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from math import isqrt,lcm
from collections import defaultdict
import json,hashlib,time,sys
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def canonical_sha(data):return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def sqrt_up(x,bits=60):
    x=F(x)
    if x<0:raise ValueError('Negative squared-norm enclosure')
    if not x:return F(0)
    scale=1<<bits;n=(x.numerator*scale*scale+x.denominator-1)//x.denominator
    v=isqrt(n)
    if v*v<n:v+=1
    return F(v,scale)
def floor_scaled(x,bits):return x.numerator*(1<<bits)//x.denominator

def apply_word(label,word):
    sign=1
    for c,p in reversed(word):
        if ((label>>p)&1)==c:return None,0
        if (label&((1<<p)-1)).bit_count()%2:sign=-sign
        label^=1<<p
    return label,sign

def labels_for(m,na,nb):
    s=m//2
    return sorted(sum(1<<(2*p) for p in a)+sum(1<<(2*p+1) for p in b)
                  for a in combinations(range(s),na) for b in combinations(range(s),nb))

def prepare_input(root=ROOT):
    t=time.monotonic();data=json.loads((root/'inputs/fixture.json').read_text());state=json.loads((root/'inputs/state.json').read_text())
    m=data['modes'];n=data['particles'];na,nb=state['spin_counts']
    if m!=16 or n!=8 or (na,nb)!=(4,4):raise ValueError('This molecular wrapper requires the declared H8 sector')
    if state['fixture_sha256']!=canonical_sha(data):raise ValueError('MPS-Hamiltonian binding')
    labels=labels_for(m,na,nb);ix={v:i for i,v in enumerate(labels)};N=len(labels)
    poly=[];den=1
    for term in data['hamiltonian']:
        word=term['word'];coeff=F(term['coefficient'])
        if len(word)>4 or any(type(c)is not int or type(p)is not int or c not in (0,1) or not 0<=p<m for c,p in word):raise ValueError('Word')
        if any(sum((2*c-1) for c,p in word if p%2==s) for s in (0,1)):raise ValueError('Hamiltonian does not preserve specified spin populations')
        poly.append((word,coeff));den=lcm(den,coeff.denominator)
    if den!=10**12:raise ValueError('Unexpected coefficient denominator')
    words=[(w,int(c*den)) for w,c in poly]
    rows=[defaultdict(int) for _ in labels]
    for j,label in enumerate(labels):
        for w,c in words:
            out,sg=apply_word(label,w)
            if sg:
                if out not in ix:raise ValueError('Unaccounted outgoing Hamiltonian action')
                rows[ix[out]][j]+=sg*c
    rows=[{j:v for j,v in row.items() if v} for row in rows]
    for i,row in enumerate(rows):
        if any(rows[j].get(i,0)!=v for j,v in row.items()):raise ValueError('Nonhermitian physical Hamiltonian')
    DD=[((s>>6)&1)+((s>>7)&1)-((s>>10)&1)-((s>>11)&1) for s in labels]
    ww=[((1,i),(0,j)) for i,j in [(6,10),(10,6),(7,11),(11,7)]]
    WR=[defaultdict(int) for _ in labels]
    for j,s in enumerate(labels):
        for w in ww:
            out,sg=apply_word(s,w)
            if sg:WR[ix[out]][j]+=sg
    WR=[dict(x) for x in WR]
    if any(abs(x)>2 for x in DD) or max(sum(abs(x) for x in row.values()) for row in WR)>2:raise ValueError('Control capacities')
    for i,row in enumerate(WR):
        if any(WR[j].get(i,0)!=v for j,v in row.items()):raise ValueError('W Hermiticity')
    # Exact MPS charge/schema validation and independent expansion.
    keys={'kind','fixture_sha256','modes','particles','spin_counts','denominator','bond_charges','tensors'}
    if set(state)!=keys or state['kind']!='integer_charge_mps_v1':raise ValueError('MPS schema')
    if state['modes']!=m or state['particles']!=n:raise ValueError('MPS sector')
    qs=state['bond_charges'];ts=state['tensors'];sd=state['denominator']
    if type(sd)is not int or sd<=0 or len(qs)!=m+1 or len(ts)!=m or qs[0]!=[[0,0]] or qs[-1]!=[[na,nb]]:raise ValueError('MPS endpoints')
    edge_tables=[]
    for site,tt in enumerate(ts):
        seen=set();rows0=defaultdict(list)
        for l,s,r,v in tt:
            if any(type(x)is not int for x in [l,s,r,v]) or s not in (0,1) or not v or not 0<=l<len(qs[site]) or not 0<=r<len(qs[site+1]) or (l,s,r)in seen:raise ValueError('MPS edge')
            seen.add((l,s,r));q=qs[site][l].copy();q[site%2]+=s
            if q!=qs[site+1][r]:raise ValueError('MPS charge flow')
            rows0[l].append((s,r,v))
        edge_tables.append(rows0)
    p=[]
    # Prefix-sharing expansion, while retaining all physical nonzero amplitudes.
    def walk(site,label,vec):
        if site==m:
            value=vec.get(0,0)
            if value:
                if label not in ix:raise ValueError('MPS outside accepted sector')
                amps[ix[label]]=value
            return
        branches=[defaultdict(int),defaultdict(int)]
        for l,a in vec.items():
            for s,r,b in edge_tables[site].get(l,()):branches[s][r]+=a*b
        for s,b in enumerate(branches):
            b={k:v for k,v in b.items() if v}
            if b:walk(site+1,label|(s<<site),b)
    amps=[0]*N;walk(0,0,{0:1});pd=sd**m
    normnum=sum(x*x for x in amps)
    if normnum<=0:raise ValueError('Zero input state')
    norm=F(normnum,pd*pd);initialD=F(sum(q*x*x for q,x in zip(DD,amps)),normnum)
    info={'sector_dimension':N,'H_stored_entries':sum(map(len,rows)),'source_terms':len(words),'MPS_expanded_nonzeros':sum(x!=0 for x in amps),'norm':str(norm),'initial_D':str(initialD),'initial_D_float':float(initialD),'initial_energy_Ha':str(F(sum(p*sum(v*amps[j] for j,v in row.items()) for p,row in zip(amps,rows)),den*normnum)),'seconds':time.monotonic()-t,'Hamiltonian_raw_sha256':sha(root/'inputs/fixture.json'),'state_raw_sha256':sha(root/'inputs/state.json')}
    return {'labels':labels,'H':rows,'D':DD,'W':WR,'p':amps,'pd':pd,'norm':norm,'den':den,'info':info}

def matvec(A,x):return [sum(a*b for a,b in zip(row,x)) for row in A]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def gram(A):
    cols=list(map(list,zip(*A)));d=len(cols);out=[[0]*d for _ in range(d)]
    for i in range(d):
        ai=cols[i]
        for j in range(i,d):out[i][j]=out[j][i]=dot(ai,cols[j])
    return out

def prepare_case(data,path):
    t=time.monotonic();c=json.loads(Path(path).read_text());N=len(data['labels'])
    if c['kind']!='fixed_number_control_trajectory_v1' or c['fixture_sha256']!=data['info']['Hamiltonian_raw_sha256'] or c['state_sha256']!=data['info']['state_raw_sha256']:raise ValueError('Case binding')
    if c['modes']!=16 or c['spin_counts']!=[4,4]:raise ValueError('Case sector')
    if c['control_operator_norm_bounds']!=['2','2']:raise ValueError('Control norm metadata')
    J=c['basis'];jd=c['basis_denominator'];hd=c['model_denominator'];zd=c['trajectory_denominator'];shift=F(c['energy_shift'])
    if shift!=F(1851,200):raise ValueError('Declared scalar phase convention')
    if len(J)!=N or type(jd)is not int or jd<=0 or not J or not J[0]:raise ValueError('Basis dimensions')
    d=len(J[0]);
    if d>128 or any(len(row)!=d or any(type(v)is not int for v in row) for row in J):raise ValueError('Basis integer entries')
    if type(hd)is not int or hd<=0 or type(zd)is not int or zd<=0:raise ValueError('Denominators')
    G=gram(J);gcap=max(F(sum(abs(x)for x in row),jd*jd) for row in G);sqrtg=sqrt_up(gcap)
    DJ=[[data['D'][i]*v for v in row] for i,row in enumerate(J)]
    DJgram=[[dot(col1,col2) for col2 in zip(*DJ)]for col1 in zip(*J)]
    # Actual full molecular action, shared across the finite control centers.
    hj=[];wj=[]
    for i in range(N):
        terms=list(data['H'][i].items());wterms=list(data['W'][i].items())
        hj.append([sum(v*J[j][k] for j,v in terms)+int(shift*data['den'])*J[i][k] for k in range(d)])
        wj.append([sum(v*J[j][k] for j,v in wterms)for k in range(d)])
    models={};gram_bits=90
    for model in c['models']:
        u=F(model['u']);v=F(model['v']);key=(str(u),str(v));h=model['h']
        if key in models or len(h)!=d or any(len(row)!=d or any(type(x)is not int for x in row)for row in h):raise ValueError('Model matrix')
        if any(h[i][j]!=h[j][i]for i in range(d)for j in range(d)):raise ValueError('Reduced Hermiticity')
        if (u*data['den']).denominator!=1 or (v*data['den']).denominator!=1:raise ValueError('Unsupported exact control grid')
        ui=int(u*data['den']);vi=int(v*data['den']);eh=[]
        for i in range(N):
            jh=[sum(J[i][l]*h[l][k] for l in range(d))for k in range(d)]
            eh.append([(hj[i][k]+ui*DJ[i][k]+vi*wj[i][k])*hd-data['den']*jh[k]for k in range(d)])
        ed=data['den']*jd*hd;eg=gram(eh)
        Dlo=[[(x*(1<<gram_bits))//(ed*ed)for x in row]for row in eg]
        models[key]={'h':h,'Dlo':Dlo,'squared_defect_denominator':ed*ed}
    return {'case':c,'J':J,'jd':jd,'hd':hd,'zd':zd,'G':G,'D':DJgram,'sqrtg':sqrtg,'models':models,'gram_bits':gram_bits,
            'dimension':d,'seconds':time.monotonic()-t,'H_applications_in_preparation':d,'input_case_sha256':sha(Path(path))}

def int_l2(x,y):return sum(a*a+b*b for a,b in zip(x,y))
def project_integral_upper(re,im,Dlo,zb,grambits):
    p=len(re)-1;d=len(re[0]);dr=[matvec(Dlo,v)for v in re];di=[matvec(Dlo,v)for v in im]
    ld=lcm(*range(1,2*p+2));ss=0
    for a in range(p+1):
        for b in range(p+1):ss+=(dot(re[a],dr[b])+dot(im[a],di[b]))*(ld//(a+b+1))
    extra=sum(abs(x)for row in re+im for x in row)**2
    answer=F(ss+ld*extra,(1<<grambits)*zb*zb*ld)
    if answer<0:raise ValueError('Inconsistent squared-defect enclosure')
    return answer

def query(data,pc,override_amp=None,override_trace=None,override_noise=None):
    t=time.monotonic();c=pc['case'];J=pc['J'];jd=pc['jd'];hd=pc['hd'];zd=pc['zd'];sg=pc['sqrtg'];N=data['norm'];p=data['p'];pd=data['pd'];d=pc['dimension']
    segments=c['segments'];err=F(0);prev=None;T=F(0);per=[];initial_error=None
    for ns,s in enumerate(segments):
        dt=F(s['duration']);key=(str(F(s['u'])),str(F(s['v'])));
        if dt<=0 or key not in pc['models']:raise ValueError('Undeclared dynamics segment')
        h=pc['models'][key]['h'];Dlo=pc['models'][key]['Dlo'];re=s['re'];im=s['im']
        if len(re)!=len(im) or not 2<=len(re)<=41 or any(len(row)!=d or any(type(x)is not int for x in row)for row in re+im):raise ValueError('Trajectory polynomial')
        pdeg=len(re)-1
        if prev is None:
            qr=matvec(J,re[0]);qi=matvec(J,im[0]);ad=jd*zd;ld=lcm(pd,ad)
            nr=[x*(ld//pd)-y*(ld//ad)for x,y in zip(p,qr)];ni=[-y*(ld//ad)for y in qi]
            initial_error=sqrt_up(F(int_l2(nr,ni),ld*ld));err+=initial_error;jump=F(0)
        else:
            dr=[x-y for x,y in zip(re[0],prev[0])];di=[x-y for x,y in zip(im[0],prev[1])]
            jump=sg*sqrt_up(F(int_l2(dr,di),zd*zd));err+=jump
        integral=project_integral_upper(re,im,Dlo,zd,pc['gram_bits']);projection=dt*sqrt_up(integral)
        # The polynomial integration error is checked independently of proposal eigensolves.
        ode=F(0)
        for l in range(pdeg+1):
            ar=matvec(h,re[l]);ai=matvec(h,im[l]);
            if l<pdeg:
                ar=[-(l+1)*hd*dt.denominator*x-dt.numerator*y for x,y in zip(im[l+1],ar)]
                ai=[ (l+1)*hd*dt.denominator*x-dt.numerator*y for x,y in zip(re[l+1],ai)]
            else:ar=[-dt.numerator*x for x in ar];ai=[-dt.numerator*x for x in ai]
            ode+=dt*sg*sqrt_up(F(int_l2(ar,ai),(zd*hd*dt.numerator)**2))/F(l+1)
        err+=projection+ode;T+=dt
        prev=([sum(row[k]for row in re)for k in range(d)],[sum(row[k]for row in im)for k in range(d)])
        per.append({'segment':ns,'end_time':str(T),'projection_defect_bound':str(projection),'integrator_defect_bound':str(ode),'jump_bound':str(jump),'accumulated_state_error':str(err)})
    re,im=prev
    def quadratic(A):return dot(re,matvec(A,re))+dot(im,matvec(A,im))
    qnorm=F(quadratic(pc['G']),(jd*zd)**2);qD=F(quadratic(pc['D']),(jd*zd)**2)
    if qnorm<=0:raise ValueError('Zero final approximation')
    value=qD/N;obs_error=2*err*(sqrt_up(N)+sqrt_up(qnorm))/N
    amp=F(c['robust_amplitude_Ha'])if override_amp is None else F(override_amp)
    tr=F(c['initial_trace_distance'])if override_trace is None else F(override_trace)
    noise=F(c.get('integrated_dephasing_rate','0')) if override_noise is None else F(override_noise)
    if amp<0 or not 0<=tr<=1 or noise<0:raise ValueError('Uncertainty bounds')
    # D,W norms <=2; |delta u|,|delta v| <=amp pointwise.
    # pure-state evolution difference <=4*T*amp; observable difference <=4 times this.
    robust_add=16*T*amp+4*tr+4*noise
    lo=value-obs_error;hi=value+obs_error
    rlo=lo-robust_add;rhi=hi+robust_add;target=F(c['target_D_upper']);targetlo=F(c.get('target_D_lower','-2'))
    return {'kind':'certified_driven_observable_v1','control_time_au':str(T),'reduced_dimension':d,
     'initial_state':'Specified rational MPS, normalized; not an assumed exact ground state',
     'initial_D':data['info']['initial_D'],'initial_D_float':data['info']['initial_D_float'],
     'state_error_bound':str(err),'state_error_float':float(err),'initial_representation_error':str(initial_error),
     'final_polynomial_norm_squared':str(qnorm),'D_center':str(value),'D_center_float':float(value),
     'D_lower':str(lo),'D_upper':str(hi),'D_interval':[float(lo),float(hi)],
     'per_control_pointwise_uncertainty_Ha':str(amp),'initial_trace_distance_radius':str(tr),
     'robust_extra_allowance':str(robust_add),'robust_D_lower':str(rlo),'robust_D_upper':str(rhi),'robust_D_interval':[float(rlo),float(rhi)],
     'target_D_lower':str(targetlo),'target_D_upper':str(target),'target_proved':bool(rhi<=target and rlo>=targetlo),'integrated_dephasing_rate_budget':str(noise),'segments':per,
     'query_seconds':time.monotonic()-t,
     'claim':'Continuous-time bound and all measurable open-loop waveform perturbations within the declared pointwise box; initial mixed states in the declared trace-distance ball within the same balanced-spin sector, and time-dependent phase-flip dephasing sum-integral within the declared rate budget. No experimentally realizable pulse, global optimum, or universal control-policy guarantee claimed.'}

def run(names,out):
    start=time.monotonic();out=Path(out)
    if out.exists():raise ValueError('Output exists')
    if ROOT.resolve() in out.resolve().parents:raise ValueError('Use an output path outside the bundle')
    out.mkdir(parents=True);data=prepare_input();(out/'input.json').write_text(json.dumps(data['info'],indent=2))
    print('inputs',data['info'],flush=True);res={}
    for name in names:
        pc=prepare_case(data,ROOT/'inputs'/f'{name}.json');print('prepared',name,pc['seconds'],flush=True)
        r=query(data,pc);r['shared_preparation_seconds']=pc['seconds'];r['input_case_sha256']=pc['input_case_sha256'];
        (out/f'{name}.json').write_text(json.dumps(r,indent=2));res[name]={k:v for k,v in r.items()if k!='segments'}
        print(name,r['D_interval'],r['robust_D_interval'],'target',r['target_proved'],'query_s',r['query_seconds'],flush=True)
    if any(x in sys.modules for x in ['numpy','scipy','pyscf','quimb','numba']):raise RuntimeError('Numerical library on accepting path')
    record={'results':res,'input':data['info'],'complete_seconds':time.monotonic()-start,'acceptance_libraries':'Python standard library','full_balanced_sector_enumerated':True,'full_fixed_N_sector_enumerated':False}
    (out/'complete.json').write_text(json.dumps(record,indent=2));return record

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--cases',nargs='+',default=['short24','long64','baseline32','short8_refusal']);p.add_argument('--out',required=True);a=p.parse_args();r=run(a.cases,a.out);print('complete',r['complete_seconds'])
