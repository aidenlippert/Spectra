"""Bounded profile/penalty discovery; exact energy replay is a separate command."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,sys,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_charged_projectors import charged_vectors
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'joint_projector'


def profiles(x,V=F(1,2)):
    a,b,p,q,d,e=map(F,x)
    return ([a,b,10-a-b,10-a-b,b,a],[p,q,5-2*p-2*q,q,p],
            [d,e,5*V-2*d-2*e,e,d])


def build(x,V=F(1,2)):
    u,t,v=profiles(x,V)
    return sector_matrices(6,F(10,3),1,u,t,V,v)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--signed-density',action='store_true')
    parser.add_argument('--telescoping',action='store_true')
    parser.add_argument('--symmetric-diagonals',type=int,default=0)
    parser.add_argument('--target-v',type=F,default=F(1,2))
    parser.add_argument('--extra-diagonal',type=int)
    args=parser.parse_args();out=OUT/'signed_density' if args.signed_density else OUT
    out.mkdir(exist_ok=True)
    if not 0<=args.symmetric_diagonals<=16:raise ValueError('At most16 diagonal corrections per probe')
    if args.symmetric_diagonals:args.telescoping=True
    if args.extra_diagonal is not None:
        if args.symmetric_diagonals!=8 or not 8<=args.extra_diagonal<120:raise ValueError('GPU-selected extra shape requires eight baseline shapes')
    Ys=[]
    if args.telescoping:
        if not args.signed_density:raise ValueError('Telescoping probe starts from signed-density checkpoint')
        if args.symmetric_diagonals:
            all_candidates=json.loads((out/'symmetry_diagonal_candidates.json').read_text())['candidates']
            candidates=all_candidates[:args.symmetric_diagonals]
            if args.extra_diagonal is not None:candidates=candidates+[all_candidates[args.extra_diagonal]]
            Ys=[{int(s):v for s,v in c['diagonal'].items()} for c in candidates]
        else:
            separator=json.loads((out/'overlap_consistency_separator.json').read_text());Y={}
            for i,j,v in separator['reflection_odd_five_site_matrix']:
                if i!=j:raise ValueError('Only diagonal separators are implemented')
                Y[i]=v
            Ys=[Y]
        out=out/(f'symmetric_diagonals_{len(Ys)}' if args.symmetric_diagonals else 'telescoping');out.mkdir(exist_ok=True)
        if args.extra_diagonal is not None:out=out/('extra_'+str(args.extra_diagonal));out.mkdir(exist_ok=True)
    if args.target_v!=F(1,2):
        if not args.signed_density:raise ValueError('Transfer probe requires signed density profiles')
        out=out/('transfer_V'+str(args.target_v).replace('/','_'));out.mkdir(exist_ok=True)
    started=time.monotonic()
    cert=json.loads((BASE/'six_site_projector/refined_certificate.json').read_text())
    if args.telescoping:cert=json.loads((OUT/'signed_density/profile_joint_r1_2_certificate.json').read_text())
    if args.extra_diagonal is not None:cert=json.loads((OUT/'signed_density/symmetric_diagonals_8/profile_joint_r1_2_certificate.json').read_text())
    if args.target_v!=F(1,2):
        cert['target']['V']=str(args.target_v);cert['local_window']['V']=str(args.target_v)
        cert['local_window']['density_profile']=[str(F(v)*2*args.target_v) for v in cert['local_window']['density_profile']]
    charged=json.loads((BASE/'charged_projector/source.json').read_text())['physical_state']
    half={int(k):v for k,v in cert['vector'].items()};hn=sum(v*v for v in half.values())
    cv,cn=charged_vectors(charged)
    local=cert['local_window'];x0=list(map(F,[local['onsite_profile'][0],local['onsite_profile'][1],local['hopping_profile'][0],local['hopping_profile'][1],local['density_profile'][0],local['density_profile'][1]]))
    original=build(x0,args.target_v);perturbed=[]
    for j in range(6):
        point=x0.copy();point[j]+=F(1,100);perturbed.append(build(point,args.target_v))
    mats=[]
    for key,(_,K,cols) in original.items():
        scale=np.sqrt([sum(a*a for a in col.values()) for col in cols])
        denom=scale[:,None]*scale[None,:];A=np.array(K,float)/denom
        deriv=np.array([(np.array(ss[key][1],float)/denom-A)*100 for ss in perturbed])
        def proj(v,n):
            w=np.array([sum(a*v.get(s,0) for s,a in col.items()) for col in cols],float)/scale/np.sqrt(float(n))
            return np.outer(w,w)
        diagonals=[]
        for Y in Ys:
            diagonal=[]
            for col in cols:
                ds={Y.get(s&1023,0)-Y.get(s>>2,0) for s in col}
                if len(ds)!=1:raise ValueError('Telescoping term breaks reflection')
                diagonal.append(ds.pop())
            diagonals.append(diagonal)
        mats.append((key,A,deriv,proj(half,hn),sum((proj(v,cn) for v in cv),np.zeros_like(A)),np.array(diagonals).reshape(len(Ys),len(cols))))
    theta_h=float(F(cert['projector_sum_ceiling'])/cert['windows']);origin=np.array(x0,float)
    rows=[]
    proposals=json.loads((OUT/'numeric_proposal.json').read_text())['rows']
    cases=[None]+proposals
    if args.telescoping:cases=[row for row in proposals if row['ratio']=='1/2']
    for row in cases:
        ratio=0 if row is None else float(F(row['ratio']));theta_j=0 if row is None else float(F(row['theta_joint']))
        calls=0;cached=None;values=None;jac=None;best=None
        def evaluate(y):
            nonlocal calls,cached,values,jac,best
            if cached is not None and np.array_equal(y,cached):return values,jac
            calls+=1
            if calls>250:raise RuntimeError('Hard250 matrix-evaluation budget exceeded')
            cached=y.copy();vals=[];grads=[]
            for key,A,deriv,P,Q,T in mats:
                matrix=A+np.einsum('i,ijk->jk',y[:6]-origin,deriv)+y[6]*P+y[7]*(P+ratio*Q)
                if args.telescoping:matrix=matrix+np.diag(y[8:-1]@T)
                ev,vec=eigh(matrix,subset_by_index=[0,0]);w=vec[:,0]
                vals.append(float(ev[0])-y[-1])
                grads.append([float(w@D@w) for D in deriv]+[float(w@P@w),float(w@(P+ratio*Q)@w)]+list(T@(w*w))+[-1.])
            values=np.array(vals);jac=np.array(grads)
            density=(min(values)+y[-1]-y[6]*theta_h-y[7]*theta_j)/5
            physical=y[0]+y[1]<=10+1e-12 and y[2]+y[3]<=2.5+1e-12
            if not args.signed_density:physical=physical and y[4]+y[5]<=1.25+1e-12
            if physical and (best is None or density>best[0]):best=(density,y.copy())
            return values,jac
        cost=np.array([0.]*6+[theta_h/5,theta_j/5]+[0.]*len(Ys)+[-.2])
        y0=np.r_[origin,float(F(cert['penalty'])),0.,float(F(cert['penalized_lower']))]
        if row is not None:
            y0[6]=float(F(row['candidate_alpha']));y0[7]=float(F(row['candidate_beta']))
            y0[8]=float(row['numerical_ell'])-1e-6
        if args.telescoping:
            seed_diagonal={int(s):F(v) for s,v in cert.get('telescoping_diagonal',{}).items()}
            seed_taus=[]
            for Y in Ys:
                ratios={seed_diagonal.get(s,F(0))/v for s,v in Y.items()}
                if len(ratios)!=1:raise ValueError('Seed correction not in selected diagonal span')
                seed_taus.append(float(ratios.pop()))
            y0=np.r_[origin,float(F(cert['penalty'])),float(F(cert['joint']['penalty'])),seed_taus,float(F(cert['penalized_lower']))]
        linear=np.zeros((3,len(y0)));linear[0,:2]=-1;linear[1,2:4]=-1;linear[2,4:6]=-1
        limits=np.array([10.,2.5,1.25])
        if args.signed_density:linear=linear[:2];limits=limits[:2]
        dbound=(-2,2) if args.signed_density else (0,1.25)
        outcome=None
        try:
            outcome=minimize(lambda y:float(cost@y),y0,jac=lambda y:cost,method='SLSQP',
                bounds=[(0,10),(0,10),(0,2.5),(0,2.5),dbound,dbound,(0,3),(0,0 if row is None else 3)]+[(-2,2)]*len(Ys)+[(-10,0)],
                constraints=[{'type':'ineq','fun':lambda y:evaluate(y)[0],'jac':lambda y:evaluate(y)[1]},
                             {'type':'ineq','fun':lambda y:limits+linear@y,'jac':lambda y:linear}],
                options={'maxiter':120,'ftol':1e-10})
        except RuntimeError as exc:
            if 'budget' not in str(exc):raise
        chosen=best[1];rational=[F(round(float(v)*10**6),10**6) for v in chosen[:8]]
        taus=[F(round(float(t)*10**6),10**6) for t in chosen[8:-1]]
        fresh=build(rational[:6],args.target_v);minimum=(float('inf'),None);max_difference=0.
        for key,A,deriv,P,Q,T in mats:
            _,K,cols=fresh[key];scale=np.sqrt([sum(a*a for a in c.values()) for c in cols]);direct=np.array(K,float)/scale[:,None]/scale[None,:]
            affine=A+np.einsum('i,ijk->jk',np.array(rational[:6],float)-origin,deriv)
            max_difference=max(max_difference,float(np.max(abs(direct-affine))))
            correction=np.array([sum(float(tau)*(Y.get(next(iter(col))&1023,0)-Y.get(next(iter(col))>>2,0)) for tau,Y in zip(taus,Ys)) for col in cols])
            if args.telescoping and np.max(abs(correction-np.array(taus,float)@T))>1e-12:raise ValueError('Fresh telescoping matrix differs')
            ell=float(eigh(direct+float(rational[6])*P+float(rational[7])*(P+ratio*Q)+np.diag(correction),eigvals_only=True,subset_by_index=[0,0])[0])
            minimum=min(minimum,(ell,key))
        if max_difference>1e-9:raise ValueError('Fresh CAR matrix disagrees with affine model')
        lower=F(math.floor(minimum[0]*10**7)-1,10**7)
        result=dict(cert,penalty=str(rational[6]),penalized_lower=str(lower),local_window=dict(local))
        onsite,hopping,density=profiles(rational[:6],args.target_v);result['local_window'].update(onsite_profile=list(map(str,onsite)),hopping_profile=list(map(str,hopping)),density_profile=list(map(str,density)))
        name='profile_half_only' if row is None else 'profile_joint_r'+row['ratio'].replace('/','_')
        if row is not None:
            result.update(kind='hubbard_projector_extension_v4',joint={'vector':charged,'windows':4,'ratio':row['ratio'],'projector_sum_ceiling':row['ceiling'],'penalty':str(rational[7])})
        if args.telescoping:
            result['kind']='hubbard_projector_extension_v5'
            correction={s:sum(tau*Y.get(s,0) for tau,Y in zip(taus,Ys)) for s in set().union(*Ys)}
            result['telescoping_diagonal']={str(s):str(v) for s,v in correction.items() if v}
            if len(result['telescoping_diagonal'])>64:result['kind']='hubbard_projector_extension_telescoping_proposal'
        bound=(lower-rational[6]*F(cert['projector_sum_ceiling'])/cert['windows']-(F(0) if row is None else rational[7]*F(row['theta_joint'])))/5
        item={'name':name,'rational_parameters':list(map(str,rational)),'local_minimum':minimum[0],'active_sector':minimum[1],'lower':str(lower),'proposed_periodic_density':str(bound),'density_float':float(bound),'matrix_evaluations':calls,'max_fresh_affine_difference':max_difference,'optimizer_success':False if outcome is None else bool(outcome.success),'optimizer_message':'hard evaluation budget' if outcome is None else str(outcome.message)}
        if args.telescoping:item['telescoping_coefficients']=list(map(str,taus))
        (out/(name+'_certificate.json')).write_text(json.dumps(result,indent=2)+'\n');rows.append(item)
        (out/'profile_proposals.json').write_text(json.dumps({'rows':rows,'seconds':time.monotonic()-started,'signed_density':args.signed_density,'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Numerical proposals only, equal250 evaluation caps, fresh exact-CAR matrix numerical check. No optimizer optimality or exact PSD claim.'},indent=2)+'\n')
        print(json.dumps(item),flush=True)


if __name__=='__main__':main()
