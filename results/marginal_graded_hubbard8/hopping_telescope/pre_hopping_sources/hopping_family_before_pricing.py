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
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value as signed_value
from experiments.marginal_quadratic_charge_telescope import local_value as quadratic_value
from experiments.marginal_charge_square_pairs import local_value as square_value
from experiments.marginal_charge_indicator_telescope import local_value as indicator_value
from experiments.marginal_range_two_density import diagonal_value


def solve_basis(matrix,rhs):
    m=len(rhs);n=len(matrix[0])
    if not 1<=n<=m<=67 or any(len(row)!=n for row in matrix):raise ValueError('Bounded rectangular dual basis required')
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
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--scale',type=float,default=.002);parser.add_argument('--determinant-anchors',action='store_true');args=parser.parse_args();out=args.directory
    if not 0<args.scale<=.2:raise ValueError('Bounded perturbation scale required')
    assert solve_basis([[1,0],[0,1],[1,1]],[2,3,5])==[2,3]
    try:solve_basis([[1],[1]],[1,2])
    except ValueError:pass
    else:raise AssertionError('Inconsistent basis accepted')
    started=time.monotonic();c=json.loads((out/'profile_joint_r1_2_certificate.json').read_text());candidates=json.loads((BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json').read_text())['candidates'];indices=[1,2,4,5,6,7,31,11,20]
    shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in indices]
    x,_,mats,half,hn,charged,cn,ratio=prepare(c,shapes)
    action=actions();hop={key:projected_matrix(cols,action) for key,a,ds,ph,q,ts,cols,*_ in mats}
    gamma=F(c['hopping_telescope'])
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
        m=float(gamma)*h+a+float(alpha)*ph+float(beta)*q+np.diag(np.array(taus)@ts)
        if eigh(m,subset_by_index=[0,0],eigvals_only=True)[0]>float(F(c['penalized_lower']))+1e-4:continue
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);exact[key]=(k,derivatives,cols,ts)
        directions=[np.zeros(64)]+[sign*np.eye(64)[j]*args.scale for j in range(64) for sign in (-1,1)]+[rng.normal(size=64)*args.scale for _ in range(16)]
        for delta in directions:
            _,vectors=eigh(m+np.einsum('i,ijk->jk',delta[:2],ds)+np.diag(delta[2:-1]@ts)+delta[-1]*h,subset_by_index=[0,0]);v=vectors[:,0]/scale;coeff=[round(float(z/max(abs(v)))*10**8) for z in v]
            vector={str(s):z*b for z,col in zip(coeff,cols) if z for s,b in col.items()};norm=sum(v*v for v in vector.values());w=np.array(coeff,float)*scale/np.sqrt(float(norm));d=np.array([w@matrix@w for matrix in ds]);row=np.r_[1,d,ts@(w*w),w@h@w,w@ph@w,w@q@w]
            if any(np.max(abs(row-old))<1e-11 for old in rows):continue
            rows.append(row);states.append((vector,norm,coeff,key));energies.append(float(w@a@w-d@np.array(x[2:4],float)))
    if len(states)>1200:raise ValueError('State-sampling budget exceeded')
    local=c['local_window'];u=list(map(F,local['onsite_profile']));density=list(map(F,local['density_profile']));range2=list(map(F,local['range_two_density_profile']))
    fixed_terms=[(quadratic_value,{k:F(v) for k,v in c['quadratic_charge_telescope'].items()}),(square_value,{k:F(v) for k,v in c['charge_square_pair_telescope'].items()}),(indicator_value,{k:F(v) for k,v in c['higher_charge_indicator_telescope'].items()})]
    def anchor(s):
        charges=[((s>>(2*i))&3).bit_count()-1 for i in range(6)]
        ph=F(half.get(s,0)**2,hn);q=ph+ratio*sum((F(v.get(s,0)**2,cn) for v in charged),F(0))
        t=[F(shape.get(s&1023,0)-shape.get(s>>2,0)) for shape in shapes]+[signed_value(s,{key:F(1)}) for key in PATTERNS]
        energy=sum((a*p*p/2 for a,p in zip(u,charges)),F(0))+sum((a*charges[i]*charges[i+1] for i,a in enumerate(density)),F(0))+diagonal_value(s,range2)+sum((fn(s,terms) for fn,terms in fixed_terms),F(0))
        return [F(1),F(0),F(0)]+t+[F(0),ph,q],energy
    if args.determinant_anchors:
        for s in range(4096):
            row,energy=anchor(s);rows.append(np.array(row,float));energies.append(float(energy));states.append(({str(s):1},1,None,None))
    if len(states)>5296:raise ValueError('Total physical candidate budget exceeded')
    arr=np.array(rows).T;rhs=[F(1)]+[F(0)]*64+[th,tj]
    result=linprog(energies,A_eq=arr[:65],b_eq=np.array(rhs[:65],float),A_ub=arr[65:],b_ub=np.array(rhs[65:],float),bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    diagnostic={'accepted':False,'scale':args.scale,'determinant_anchors':args.determinant_anchors,'states':len(states),'active_sectors':list(exact),'status':result.message,'seconds':time.monotonic()-started}
    if result.success:
        values=list(result.x)+list(np.array(rhs[65:],float)-arr[65:]@result.x);support=[i for i,v in enumerate(values) if v>1e-10];diagnostic.update(numerical_family_upper=float(result.fun)/5,basis_size=len(support))
        erows=[];eenergies=[]
        for index in support:
            if index>=len(states):erows.append([F(int(j==65+index-len(states))) for j in range(67)]);eenergies.append(F(0));continue
            vector,norm,coeff,key=states[index]
            if key is None:
                row,energy=anchor(int(next(iter(vector))));erows.append(row);eenergies.append(energy);continue
            k,derivatives,cols,ts=exact[key]
            def expectation(matrix):return sum((F(a*b)*matrix[i][j] for i,a in enumerate(coeff) if a for j,b in enumerate(coeff) if b and matrix[i][j]),F(0))/norm
            d=[expectation(matrix) for matrix in derivatives];t=[F(sum(int(ts[j,i])*z*z*sum(b*b for b in cols[i].values()) for i,z in enumerate(coeff)),norm) for j in range(61)]
            def fidelity(source,snorm):return F(sum(v*source.get(int(s),0) for s,v in vector.items())**2,norm*snorm)
            ph=fidelity(half,hn);q=ph+ratio*sum((fidelity(v,cn) for v in charged),F(0));erows.append([F(1)]+d+t+[expectation(hop[key]),ph,q]);eenergies.append(expectation(k)-sum(a*b for a,b in zip(x[2:4],d)))
        try:
            weights=solve_basis([[row[i] for row in erows] for i in range(67)],rhs)
            if min(weights)<0:raise ValueError('Exact basis has negative weights')
            upper=sum(w*e for w,e in zip(weights,eenergies))/5
            proposal={'kind':'joint_hopping_telescope_family_proposal_v1','half_vector':c['vector'],'charged_vector':c['joint']['vector'],'ratio':str(ratio),'theta_half':str(th),'theta_joint':str(tj),'diagonal_shapes':[{str(s):str(v) for s,v in shape.items()} for shape in shapes],'mixture':[{'weight':str(w),'vector':states[i][0]} for i,w in zip(support,weights) if i<len(states) and w],'W':c['target']['W'],'range_two_density_profile':c['local_window']['range_two_density_profile'],'proposed_periodic_family_upper':str(upper),'scope':'Untrusted exact rectangular-basis proposal; every constraint and energy expectation requires production physical replay.'}
            (out/'diagonal_family_limit_proposal.json').write_text(json.dumps(proposal,indent=2)+'\n');diagnostic.update(rational_family_upper=str(upper),rational_family_upper_float=float(upper),mixture_sources=len(proposal['mixture']))
        except ValueError as exc:diagnostic['exact_basis_failure']=str(exc)
    diagnostic['seconds']=time.monotonic()-started;(out/'diagonal_family_limit_diagnostic.json').write_text(json.dumps(diagnostic,indent=2)+'\n');print(json.dumps(diagnostic))


if __name__=='__main__':main()
