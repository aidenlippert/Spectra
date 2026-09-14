"""Nonaccepting bounded dual search for the complete signed-charge chart."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from signed_charge_numeric import prepare,BASE
from experiments.marginal_hopping_telescope import actions,projected_matrix
from experiments.marginal_spin_telescope import LABELS,actions as spin_actions
from experiments.marginal_spectator_hopping import LABELS as SPECTATORS,actions as spectator_actions
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value as signed_value
from experiments.marginal_quadratic_charge_telescope import local_value as quadratic_value
from experiments.marginal_charge_square_pairs import local_value as square_value
from experiments.marginal_charge_indicator_telescope import local_value as indicator_value
from experiments.marginal_range_two_density import diagonal_value


def solve_basis(matrix,rhs):
    m=len(rhs);n=len(matrix[0])
    if not 1<=n<=m<=85 or any(len(row)!=n for row in matrix):raise ValueError('Bounded rectangular dual basis required')
    a=[list(map(F,row))+[F(b)] for row,b in zip(matrix,rhs)]
    for j in range(n):
        pivot=next((i for i in range(j,m) if a[i][j]),None)
        if pivot is None:raise ValueError('Dependent proposed basis columns')
        a[j],a[pivot]=a[pivot],a[j];v=a[j][j];a[j]=[z/v for z in a[j]]
        for i in range(m):
            if i!=j and a[i][j]:
                v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
    if any(row[-1] for row in a[n:]):raise ValueError('Exact proposed basis is inconsistent')
    return [a[j][-1] for j in range(n)]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--extra-atoms',type=Path,required=True);parser.add_argument('--scale',type=float,default=.002);parser.add_argument('--determinant-anchors',action='store_true');parser.add_argument('--pricing-rounds',type=int,default=0);parser.add_argument('--pricing-radius',type=float,default=0);parser.add_argument('--pricing-eigenvectors',type=int,choices=(1,2),default=1);args=parser.parse_args();out=args.directory
    if not 0<=args.pricing_radius<=.2:raise ValueError('Pricing radius must be zero or at most .2')
    if args.pricing_eigenvectors==2 and args.pricing_rounds>40:raise ValueError('Two-eigenvector pricing cap40 preserves the total candidate budget')
    if not 0<=args.pricing_rounds<=80:raise ValueError('Pricing round cap80 exceeded')
    if not 0<args.scale<=.2:raise ValueError('Bounded perturbation scale required')
    assert solve_basis([[1,0],[0,1],[1,1]],[2,3,5])==[2,3]
    try:solve_basis([[1],[1]],[1,2])
    except ValueError:pass
    else:raise AssertionError('Inconsistent basis accepted')
    started=time.monotonic();c=json.loads((out/'profile_joint_r1_2_certificate.json').read_text());candidates=json.loads((BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json').read_text())['candidates'];indices=[1,2,4,5,6,7,31,11,20]
    shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in indices]
    x,_,mats,half,hn,charged,cn,ratio=prepare(c,shapes)
    action=[actions()]+[spin_actions({label:1}) for label in LABELS]+[spectator_actions({label:1}) for label in SPECTATORS];hop={key:[projected_matrix(cols,item) for item in action] for key,a,ds,ph,q,ts,cols,*_ in mats}
    gamma=[F(c['hopping_telescope'])]+[F(c['spin_telescope'].get(label,0)) for label in LABELS]+[F(c['spectator_hopping'].get(label,0)) for label in SPECTATORS]
    diagonal={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};taus=[]
    for shape in shapes:
        values={diagonal.get(s,F(0))/v for s,v in shape.items()}
        if len(values)!=1:raise ValueError('Fixed sparse span mismatch')
        taus.append(float(values.pop()))
    taus+= [float(F(c['signed_charge_telescope'].get(key,0))) for key in PATTERNS]
    th=F(c['projector_sum_ceiling'])/c['windows'];tj=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows'];alpha=F(c['penalty']);beta=F(c['joint']['penalty'])
    states=[];rows=[];energies=[];exact={};rng=np.random.default_rng(20260912)
    for key,a,ds,ph,q,ts,cols,k,derivatives in mats:
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);h=np.array(hop[key],float)/scale[:,None]/scale[None,:]
        m=np.einsum('i,ijk->jk',np.array(gamma,float),h)+a+float(alpha)*ph+float(beta)*q+np.diag(np.array(taus)@ts)
        if eigh(m,subset_by_index=[0,0],eigvals_only=True)[0]>float(F(c['penalized_lower']))+1e-4:continue
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);exact[key]=(k,derivatives,cols,ts)
        directions=[np.zeros(82)]+[sign*np.eye(82)[j]*args.scale for j in range(82) for sign in (-1,1)]+[rng.normal(size=82)*args.scale for _ in range(16)]
        for delta in directions:
            _,vectors=eigh(m+np.einsum('i,ijk->jk',delta[:2],ds)+np.diag(delta[2:63]@ts)+np.einsum('i,ijk->jk',delta[63:82],h),subset_by_index=[0,0]);v=vectors[:,0]/scale;coeff=[round(float(z/max(abs(v)))*10**8) for z in v]
            vector={str(s):z*b for z,col in zip(coeff,cols) if z for s,b in col.items()};norm=sum(v*v for v in vector.values());w=np.array(coeff,float)*scale/np.sqrt(float(norm));d=np.array([w@matrix@w for matrix in ds]);row=np.r_[1,d,ts@(w*w),[w@matrix@w for matrix in h],w@ph@w,w@q@w]
            if any(np.max(abs(row-old))<1e-11 for old in rows):continue
            rows.append(row);states.append((vector,norm,coeff,key));energies.append(float(w@a@w-d@np.array(x[2:4],float)))
    # Add rounded physical eigenvectors from the bounded conic dual, plus
    # small deterministic tangent perturbations. Nothing is accepted here.
    atom_data=json.loads(args.extra_atoms.read_text())
    atoms=sorted(atom_data['atoms'],key=lambda item:item['dual_weight'],reverse=True)
    atoms=[item for item in atoms[:16] if item['dual_weight']>1e-6]
    by_key={item[0]:item for item in mats}
    added_atoms=0;conic_rows=[];conic_energies=[]
    for atom in atoms:
        key=tuple(atom['sector']);key,a,ds,ph,q,ts,cols,k,derivatives=by_key[key]
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols])
        h=np.array(hop[key],float)/scale[:,None]/scale[None,:]
        raw=np.array(atom['coordinates']);assert len(raw)==len(cols)
        directions=[raw]+[raw+epsilon*rng.normal(size=len(raw))/np.sqrt(len(raw)) for epsilon in (1e-4,1e-3) for _ in range(12)]
        exact[key]=(k,derivatives,cols,ts)
        for direction_index,raw_vector in enumerate(directions):
            v=raw_vector/scale;coeff=[round(float(z/max(abs(v)))*10**8) for z in v]
            vector={str(s):z*b for z,col in zip(coeff,cols) if z for s,b in col.items()}
            norm=sum(v*v for v in vector.values());w=np.array(coeff,float)*scale/np.sqrt(float(norm))
            d=np.array([w@matrix@w for matrix in ds]);row=np.r_[1,d,ts@(w*w),[w@matrix@w for matrix in h],w@ph@w,w@q@w]
            rows.append(row);states.append((vector,norm,coeff,key));energies.append(float(w@a@w-d@np.array(x[2:4],float)));added_atoms+=1
            if direction_index==0:conic_rows.append(5*atom["dual_weight"]*row);conic_energies.append(5*atom["dual_weight"]*energies[-1])
    print(json.dumps({'conic_atoms':len(atoms),'extra_sampled_vectors':added_atoms,'sampled_states':len(states)}),flush=True)
    moment=sum(conic_rows);residual=moment-np.r_[1,np.zeros(82),float(th),float(tj)]
    result={'accepted':False,'trace':float(moment[0]),'moments':moment.tolist(),'residuals':residual.tolist(),'approximate_family_upper':sum(conic_energies)/5,'largest_residuals':sorted(enumerate(residual),key=lambda x:abs(x[1]),reverse=True)[:12]}
    (out/'conic_atom_moment_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('moments','residuals')}));return

    local=c['local_window'];u=list(map(F,local['onsite_profile']));density=list(map(F,local['density_profile']));range2=list(map(F,local['range_two_density_profile']))
    fixed_terms=[(quadratic_value,{k:F(v) for k,v in c['quadratic_charge_telescope'].items()}),(square_value,{k:F(v) for k,v in c['charge_square_pair_telescope'].items()}),(indicator_value,{k:F(v) for k,v in c['higher_charge_indicator_telescope'].items()})]
    def anchor(s):
        charges=[((s>>(2*i))&3).bit_count()-1 for i in range(6)]
        ph=F(half.get(s,0)**2,hn);q=ph+ratio*sum((F(v.get(s,0)**2,cn) for v in charged),F(0))
        t=[F(shape.get(s&1023,0)-shape.get(s>>2,0)) for shape in shapes]+[signed_value(s,{key:F(1)}) for key in PATTERNS]
        energy=sum((a*p*p/2 for a,p in zip(u,charges)),F(0))+sum((a*charges[i]*charges[i+1] for i,a in enumerate(density)),F(0))+diagonal_value(s,range2)+sum((fn(s,terms) for fn,terms in fixed_terms),F(0))
        return [F(1),F(0),F(0)]+t+[F(item[s].get(s,0)) for item in action]+[ph,q],energy
    if args.determinant_anchors:
        for s in range(4096):
            row,energy=anchor(s);rows.append(np.array(row,float));energies.append(float(energy));states.append(({str(s):1},1,None,None))
    if len(states)>5296:raise ValueError('Total physical candidate budget exceeded')
    arr=np.array(rows).T;rhs=[F(1)]+[F(0)]*82+[th,tj]
    result=linprog(energies,A_eq=arr[:83],b_eq=np.array(rhs[:83],float),A_ub=arr[83:],b_ub=np.array(rhs[83:],float),bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
    pricing=[]
    seed_dual=np.r_[float(F(c['penalized_lower'])),-np.array(x[2:4],float),-np.array(taus),-np.array(gamma,float),-float(alpha),-float(beta)]
    for iteration in range(args.pricing_rounds):
        if not result.success:break
        dual=np.r_[result.eqlin.marginals,result.ineqlin.marginals];added=0;lowest=float('inf')
        if args.pricing_radius:dual=seed_dual+np.clip(dual-seed_dual,-args.pricing_radius,args.pricing_radius)
        for key,a,ds,ph,q,ts,cols,k,derivatives in mats:
            scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);h=np.array(hop[key],float)/scale[:,None]/scale[None,:]
            price=a-np.einsum('i,ijk->jk',np.array(x[2:4],float)+dual[1:3],ds)-np.diag(dual[3:64]@ts)-np.einsum('i,ijk->jk',dual[64:83],h)-dual[83]*ph-dual[84]*q
            ev,vectors=eigh(price,subset_by_index=[0,min(args.pricing_eigenvectors,len(price))-1]);lowest=min(lowest,float(ev[0])-dual[0])
            for column,value in enumerate(ev):
                if float(value)-dual[0]>=-1e-9:continue
                v=vectors[:,column]/scale;coeff=[round(float(z/max(abs(v)))*10**8) for z in v]
                vector={str(s):z*b for z,col in zip(coeff,cols) if z for s,b in col.items()};norm=sum(v*v for v in vector.values());w=np.array(coeff,float)*scale/np.sqrt(float(norm));d=np.array([w@matrix@w for matrix in ds]);row=np.r_[1,d,ts@(w*w),[w@matrix@w for matrix in h],w@ph@w,w@q@w]
                exact[key]=(k,derivatives,cols,ts);rows.append(row);states.append((vector,norm,coeff,key));energies.append(float(w@a@w-d@np.array(x[2:4],float)));added+=1
        pricing.append({'round':iteration,'sampled_upper':float(result.fun)/5,'minimum_reduced_eigenvalue':lowest,'added':added})
        if not added:break
        if len(states)>13000:raise ValueError('Pricing physical candidate cap exceeded')
        arr=np.array(rows).T
        result=linprog(energies,A_eq=arr[:83],b_eq=np.array(rhs[:83],float),A_ub=arr[83:],b_ub=np.array(rhs[83:],float),bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-10,'primal_feasibility_tolerance':1e-10})
    diagnostic={'accepted':False,'proposal_written':False,'lp_tolerance':1e-10,'support_threshold':1e-14,'scale':args.scale,'extra_atom_file':str(args.extra_atoms),'extra_sampled_vectors':added_atoms,'pricing':pricing,'pricing_eigenvectors':args.pricing_eigenvectors,'pricing_radius':args.pricing_radius,'determinant_anchors':args.determinant_anchors,'states':len(states),'active_sectors':list(exact),'status':result.message,'seconds':time.monotonic()-started}
    if result.success:
        values=list(result.x)+list(np.array(rhs[83:],float)-arr[83:]@result.x);support=[i for i,v in enumerate(values) if v>1e-14];diagnostic.update(numerical_family_upper=float(result.fun)/5,basis_size=len(support))
        erows=[];eenergies=[]
        for index in support:
            if index>=len(states):erows.append([F(int(j==83+index-len(states))) for j in range(85)]);eenergies.append(F(0));continue
            vector,norm,coeff,key=states[index]
            if key is None:
                row,energy=anchor(int(next(iter(vector))));erows.append(row);eenergies.append(energy);continue
            k,derivatives,cols,ts=exact[key]
            def expectation(matrix):return sum((F(a*b)*matrix[i][j] for i,a in enumerate(coeff) if a for j,b in enumerate(coeff) if b and matrix[i][j]),F(0))/norm
            d=[expectation(matrix) for matrix in derivatives];t=[F(sum(int(ts[j,i])*z*z*sum(b*b for b in cols[i].values()) for i,z in enumerate(coeff)),norm) for j in range(61)]
            def fidelity(source,snorm):return F(sum(v*source.get(int(s),0) for s,v in vector.items())**2,norm*snorm)
            ph=fidelity(half,hn);q=ph+ratio*sum((fidelity(v,cn) for v in charged),F(0));erows.append([F(1)]+d+t+[expectation(matrix) for matrix in hop[key]]+[ph,q]);eenergies.append(expectation(k)-sum(a*b for a,b in zip(x[2:4],d)))
        try:
            weights=solve_basis([[row[i] for row in erows] for i in range(85)],rhs)
            if min(weights)<0:raise ValueError('Exact basis has negative weights')
            upper=sum(w*e for w,e in zip(weights,eenergies))/5
            proposal={'kind':'joint_spectator_hopping_family_proposal_v1','half_vector':c['vector'],'charged_vector':c['joint']['vector'],'ratio':str(ratio),'theta_half':str(th),'theta_joint':str(tj),'diagonal_shapes':[{str(s):str(v) for s,v in shape.items()} for shape in shapes],'mixture':[{'weight':str(w),'vector':states[i][0]} for i,w in zip(support,weights) if i<len(states) and w],'W':c['target']['W'],'range_two_density_profile':c['local_window']['range_two_density_profile'],'proposed_periodic_family_upper':str(upper),'scope':'Untrusted exact rectangular-basis proposal; every constraint and energy expectation requires production physical replay.'}
            (out/'diagonal_family_limit_proposal.json').write_text(json.dumps(proposal,indent=2)+'\n');diagnostic.update(proposal_written=True,rational_family_upper=str(upper),rational_family_upper_float=float(upper),mixture_sources=len(proposal['mixture']))
        except ValueError as exc:diagnostic['exact_basis_failure']=str(exc)
    diagnostic['seconds']=time.monotonic()-started;(out/'diagonal_family_limit_diagnostic.json').write_text(json.dumps(diagnostic,indent=2)+'\n');print(json.dumps(diagnostic))


if __name__=='__main__':main()
