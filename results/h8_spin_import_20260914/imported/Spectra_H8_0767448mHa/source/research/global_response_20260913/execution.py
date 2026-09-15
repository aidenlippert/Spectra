"""Matched physical-vector execution; all determinant use is diagnostic."""
from fractions import Fraction as F
import json
import resource
import sys
import time

from research.compact_response_20260913 import program,closure
from research.composable_response_20260913 import joint
from research.global_response_20260913 import global_program as gp


class Counter:
    def __init__(self,matrix):
        self.matrix=matrix;self.H_rhs=0;self.projector_calls=0;self.recurrence_vector_entries=0
    def H(self,v):
        if v.ndim!=1:raise ValueError('Every RHS column must be counted separately')
        self.H_rhs+=1;return self.matrix@v
    def project(self,mask,v):
        self.projector_calls+=1;return mask*v


def response(action,v,r,counter,order=None):
    import numpy as np
    k=r['order'] if order is None else order
    delta=float(F(r['delta_Ha']));M=float(F(r['M_Ha']));c=(M+delta)/2;a=(M-delta)/2;z=c/a
    previous=np.zeros_like(v);current=v/a;tp=1.;tc=z
    for _ in range(1,k):
        Z=(c*current-action(current))/a
        previous,current=current,2*Z-previous+2*tc/a*v
        tp,tc=tc,2*z*tc-tp;counter.recurrence_vector_entries+=len(v)
    return current/tc


def retained(action,P,Q,v,r,counter,factored=False):
    Hv=action(counter.project(P,v));W=counter.project(Q,Hv)
    def D(x):return counter.project(Q,action(counter.project(Q,x)))
    if factored:
        # F_k = [T_2k(z0)/(T_2k(z0)+1)] p_2k. Noncommutative
        # order is unchanged because all powers here are of this same D.
        d=F(r['delta_Ha']);M=F(r['M_Ha']);T=program.chebyshev(2*r['order'],(M+d)/(M-d))
        FW=float(T/(T+1))*response(D,W,r,counter,2*r['order'])
    else:
        pW=response(D,W,r,counter);DpW=D(pW);pDpW=response(D,DpW,r,counter)
        FW=2*pW-pDpW
    return counter.project(P,Hv)-counter.project(P,action(FW))


def dense_retained(M,P,Q,r):
    import numpy as np
    A=M[np.ix_(P,P)]
    if not Q:return A,{'dimension':0}
    D=M[np.ix_(Q,Q)];B=M[np.ix_(P,Q)];e,V=np.linalg.eigh(D)
    d=float(F(r['delta_Ha']));maximum=float(F(r['M_Ha']));c=(maximum+d)/2;a=(maximum-d)/2
    if e[0]<d-1e-9 or e[-1]>maximum+1e-9:raise ValueError('Physical diagnostic spectrum violates certified response interval')
    coeff=np.zeros(r['order']+1);coeff[-1]=1
    residual=np.polynomial.chebyshev.chebval((c-e)/a,coeff)/np.polynomial.chebyshev.chebval(c/a,coeff)
    fd=(1-residual*residual)/e
    K=A-(B@V*fd)@(B@V).T
    return K,{'dimension':len(Q),'numeric_min_eigenvalue_Ha':float(e[0]),'numeric_max_eigenvalue_Ha':float(e[-1])}


def run():
    import numpy as np
    start=time.monotonic();data,tail,ref=program.load_case('h6');groups,den,preparation=closure.blocks(data)
    first=json.loads((program.OUT/'h6_program.json').read_text());old=json.loads((joint.OUT/'h6_joint_response.json').read_text())
    second=old['second_response'];new=json.loads((gp.OUT/'h6_global.json').read_text());rnew=new['response']
    if first['target_Ha']!=new['target_Ha']:raise ValueError('Execution targets must match exactly')
    eta1=float(F(program.check(data,tail,first)['exact_program_residual_penalty_Ha']))
    eta2=float(F(second['residual_penalty_Ha']));etanew=float(F(rnew['residual_penalty_Ha']));b=float(F(first['target_Ha']))
    if etanew>eta1+eta2:raise ValueError('New response spends a larger total allowance')
    diagnostics=[];chosen=None
    for group in groups:
        states=group['states'];M=np.array(group['H'],dtype=float)/float(den)-b*np.eye(len(states))
        Q1=[i for i,s in enumerate(states) if ((s>>10)&3)==3];P1=[i for i in range(len(states)) if i not in Q1]
        R2=[i for i in P1 if ((states[i]>>8)&3)==3];P2=[i for i in P1 if i not in R2]
        K1,d1=dense_retained(M,P1,Q1,first);L1=K1-eta1*np.eye(len(P1))
        ip=[P1.index(i) for i in P2];iq=[P1.index(i) for i in R2]
        K2,d2=dense_retained(L1,ip,iq,second);terminal=K2-eta2*np.eye(len(P2))
        Kg,dg=dense_retained(M,P2,Q1+R2,rnew);tg=Kg-etanew*np.eye(len(P2))
        e=np.linalg.eigvalsh(terminal);eg=np.linalg.eigvalsh(tg)
        diagnostics.append({'block_key':group['key'],'terminal_dimension':len(P2),
            'old_terminal_min_float_Ha':float(e[0]),'global_terminal_min_float_Ha':float(eg[0]),
            'old_terminal_eigenvalues_below_0_02_Ha':int(np.sum(e<.02)),'lowest_old_terminal_eigenvalues_Ha':e[:5].tolist(),
            'D1_spectrum':d1,'D2_spectrum':d2,'combined_D_spectrum':dg,'terminal_matrix_entries_constructed':2*len(P2)**2})
        if 63 in states:chosen=(states,M,Q1,R2,P2,terminal,tg)
    states,M,Q1,R2,P2,terminal,tg=chosen;n=len(states);v=np.zeros(n);v[states.index(63)]=1
    q1=np.array([i in Q1 for i in range(n)],dtype=float);p1=1-q1
    q2=np.array([i in R2 for i in range(n)],dtype=float);p2=p1-q2
    rows=[]
    for kind in ('original_nested','equivalent_factored','simultaneous_global'):
        counter=Counter(M);fused=kind=='equivalent_factored'
        def L1(x):return retained(counter.H,p1,q1,x,first,counter,fused)-eta1*counter.project(p1,x)
        t=time.monotonic()
        if kind=='simultaneous_global':out=retained(counter.H,p2,q1+q2,v,rnew,counter,True)-etanew*counter.project(p2,v);oracle=tg
        else:out=retained(L1,p2,q2,v,second,counter,fused)-eta2*counter.project(p2,v);oracle=terminal
        elapsed=time.monotonic()-t;expected=np.zeros(n);expected[P2]=oracle@v[P2]
        error=float(np.max(np.abs(out-expected)))
        if error>1e-9:raise ValueError('Physical action disagrees with independent spectral diagnostic')
        expected_calls=2*rnew['order']+1 if kind=='simultaneous_global' else (2*first['order']+1)*(2*second['order']+1)
        if counter.H_rhs!=expected_calls:raise ValueError('Primitive action accounting mismatch')
        rows.append({'method':kind,'seconds':elapsed,'base_H_RHS_applications':counter.H_rhs,
            'projector_calls':counter.projector_calls,'recurrence_vector_entries':counter.recurrence_vector_entries,
            'dense_H_scalar_multiply_add_pairs':counter.H_rhs*n*n,'max_abs_error_vs_spectral_float_oracle':error,
            'vector_dimension':n,'coefficient_precision':'float64; accepting certificates use exact rationals',
            'dense_many_body_H_entries':n*n,'many_body_vector_support_representation':'full conserved-block vector, diagnostic only'})
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
    result={'target_Ha':first['target_Ha'],'first_order':first['order'],'second_order':second['order'],'global_order':rnew['order'],
        'nested_total_uniform_allowance_Ha':eta1+eta2,'global_uniform_allowance_Ha':etanew,
        'preparation':preparation,'all_terminal_block_diagnostics':diagnostics,'applications':rows,
        'action_count_reduction':rows[0]['base_H_RHS_applications']/rows[2]['base_H_RHS_applications'],
        'wall_time_ratio_nested_over_global':rows[0]['seconds']/rows[2]['seconds'],'peak_RSS_bytes':peak,
        'wall_seconds':time.monotonic()-start,'scope':'Charged enumerated physical diagnostic and matched float64 execution, not determinant-free physical-vector execution. Spectral floating values are diagnostics; terminal positivity is independently reference-certified.'}
    (gp.OUT/'execution.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'applications':rows,'action_count_reduction':result['action_count_reduction'],'seconds':result['wall_seconds']},indent=2))


if __name__=='__main__':run()
