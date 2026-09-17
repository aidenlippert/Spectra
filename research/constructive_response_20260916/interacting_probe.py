"""Charged, enumerated Hubbard diagnostics; never an exact accepting checker.

Compare bare and locally dressed retained projectors of the SAME Hamiltonian.
Angles use only U and the isolated rung. Eigenvectors and sector matrices are
explicit numerical diagnostics, not a claimed compact construction.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import acosh, ceil, log, sqrt
from pathlib import Path
import time

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh


RATIONAL_HALF_ANGLE = F(116434, 10**6)


def balanced_basis(sites):
    if sites < 2 or sites % 2:
        raise ValueError('Even site count required')
    choices = list(combinations(range(sites), sites // 2))
    return sorted(sum(1 << (2*i) for i in up) + sum(1 << (2*i+1) for i in down)
                  for up in choices for down in choices)


def hop(label, target, source):
    """Exact c_target^dagger c_source action in ascending interleaved order."""
    if not label & (1 << source) or label & (1 << target):
        return None
    sign = -1 if (label & ((1 << source)-1)).bit_count() % 2 else 1
    intermediate = label ^ (1 << source)
    if (intermediate & ((1 << target)-1)).bit_count() % 2:
        sign = -sign
    return intermediate | (1 << target), sign


def doublons(label, sites):
    return sum(((label >> (2*i)) & 3) == 3 for i in range(sites))


def ladder_edges(rungs, coupling):
    lam = F(coupling)
    if not 0 <= lam <= 1:
        raise ValueError('Continuation is in [0,1]')
    result = [(2*i, 2*i+1, F(1)) for i in range(rungs)]
    result += [(2*i+s, 2*(i+1)+s, lam) for i in range(rungs-1) for s in (0,1)]
    return [(i,j,t) for i,j,t in result if t]


def rational_rows(rungs, coupling, U=F(8)):
    sites = 2*rungs
    labels = balanced_basis(sites)
    lookup = {v:i for i,v in enumerate(labels)}
    rows = [dict() for _ in labels]
    edges = ladder_edges(rungs, coupling)
    for col, label in enumerate(labels):
        value = U * doublons(label, sites)
        if value:
            rows[col][col] = value
        for a,b,t in edges:
            for s in (0,1):
                for target,source in ((2*a+s,2*b+s),(2*b+s,2*a+s)):
                    moved = hop(label,target,source)
                    if moved is not None:
                        out, sign = moved
                        row = lookup[out]
                        rows[row][col] = rows[row].get(col,F(0)) - t*sign
    for i,row in enumerate(rows):
        for j,value in row.items():
            if rows[j].get(i,F(0)) != value:
                raise AssertionError('Non-Hermitian exact Hamiltonian')
    return labels,rows


def numerical(rows):
    rr=[]; cc=[]; vv=[]
    for i,row in enumerate(rows):
        for j,x in row.items():
            if x:
                rr.append(i); cc.append(j); vv.append(float(x))
    return coo_matrix((vv,(rr,cc)),shape=(len(rows),len(rows))).tocsr()


def local_rotation(r=RATIONAL_HALF_ANGLE):
    """Exact 16-state gate. U|S> = c|S> + s|D+>, preserving N and Sz."""
    r=F(r); c=(1-r*r)/(1+r*r); s=2*r/(1+r*r)
    # Twice the projectors/cross products avoids introducing sqrt(2).
    singlet={9:1,6:-1}; ionic={3:1,12:1}
    g=[[F(i==j) for j in range(16)] for i in range(16)]
    for vector in (singlet,ionic):
        for i,x in vector.items():
            for j,y in vector.items():
                g[i][j] += (c-1)*F(x*y,2)
    for i,x in ionic.items():
        for j,y in singlet.items():
            g[i][j] += s*F(x*y,2)
            g[j][i] -= s*F(x*y,2)
    for i in range(16):
        for j in range(16):
            if sum(g[k][i]*g[k][j] for k in range(16)) != int(i==j):
                raise AssertionError('Rational local gate is not orthogonal')
    return g


def lifted_gate(labels, rung, gate):
    lookup={label:i for i,label in enumerate(labels)}
    rr=[]; cc=[]; vv=[]; shift=4*rung; mask=15 << shift
    for j,label in enumerate(labels):
        local=(label >> shift)&15
        for target in range(16):
            x=gate[target][local]
            if x:
                out=(label & ~mask) | (target << shift)
                rr.append(lookup[out]); cc.append(j); vv.append(float(x))
    return coo_matrix((vv,(rr,cc)),shape=(len(labels),len(labels))).tocsr()


def smallest(h, vector=False):
    if h.shape[0] < 20:
        values,vectors=np.linalg.eigh(h.toarray())
        return (float(values[0]),vectors[:,0]) if vector else float(values[0])
    # Deterministic non-symmetry-restricted start. Solver output is untrusted.
    v0=np.sin(np.arange(h.shape[0],dtype=float)+0.731)
    values,vectors=eigsh(h,k=1,which='SA',v0=v0,tol=2e-12,maxiter=20000)
    return (float(values[0]),vectors[:,0]) if vector else float(values[0])


def largest(h):
    if h.shape[0] < 20:
        return float(np.linalg.eigvalsh(h.toarray())[-1])
    v0=np.cos(np.arange(h.shape[0],dtype=float)+0.19)
    return float(eigsh(h,k=1,which='LA',v0=v0,tol=2e-12,
                       maxiter=20000,return_eigenvectors=False)[0])


def opnorm_small_right(x):
    gram=np.asarray(x.T@x)
    return sqrt(max(0,float(np.linalg.eigvalsh((gram+gram.T)/2)[-1])))


def analyze(h, pidx, psi, e0, target_margin, max_steps):
    start=time.monotonic()
    n=h.shape[0]; keep=np.zeros(n,dtype=bool);keep[pidx]=True
    qidx=np.flatnonzero(~keep)
    c=h[qidx,:][:,qidx].tocsr()
    b=h[qidx,:][:,pidx].toarray()
    a=h[pidx,:][:,pidx].toarray()
    p=float(psi[pidx]@psi[pidx]);q=1-p
    cmin=smallest(c);cmax=largest(c)
    energy=e0-target_margin
    A=c-energy*diags(np.ones(len(qidx)),format='csr')
    delta=cmin-energy;upper=cmax-energy
    if delta <= 0:
        raise RuntimeError('Numerical eliminated block is not positive')
    u=psi[pidx]; numerator=float(u @ (a-e0*np.eye(len(pidx))) @ u)
    ground_residual=float(np.linalg.norm(h@psi-e0*psi))
    gap_upper=numerator/q if q > 1e-12 else None
    if gap_upper is not None and cmin-e0 > gap_upper+1e-7:
        raise AssertionError('Retained-weight Rayleigh bound violated')
    response_norm_floor=sqrt(max(0,q/p))
    g=opnorm_small_right(b)
    condition=upper/delta
    # Scalar interval degree is only a floating diagnostic (not accepted gap).
    if upper-delta < 1e-13:
        chebyshev_degree=1
    else:
        need=g/sqrt(delta*target_margin)
        chebyshev_degree=max(1,ceil(acosh(max(1,need))/acosh((upper+delta)/(upper-delta))))
    diagonal=A.diagonal()
    if np.min(diagonal) <= 0:
        raise AssertionError('Positive matrix requires positive diagonal')
    off=A-diags(diagonal,format='csr')
    invsqrt=diags(1/np.sqrt(diagonal),format='csr')
    relative=invsqrt@off@invsqrt
    rho=max(abs(smallest(relative)),abs(largest(relative)))
    X=np.zeros_like(b)
    sample_steps=set([0,1,2,4,8,16,32,64,max_steps])
    history=[]; first=None
    for step in range(max_steps+1):
        if step in sample_steps:
            residual=b-A@X
            residual_norm=opnorm_small_right(residual)
            allowance=residual_norm**2/delta
            history.append({'iterations':step,'residual_operator_norm':residual_norm,
                            'gap_based_response_error_bound':allowance})
            if first is None and allowance <= target_margin:
                first=step
        if step < max_steps:
            X=(b-off@X)/diagonal[:,None]
            if not np.all(np.isfinite(X)):
                break
    return {
        'retained_dimension':len(pidx),'eliminated_dimension':len(qidx),
        'retained_ground_weight':p,'ground_residual_l2':ground_residual,
        'eliminated_ground_energy':cmin,'eliminated_gap_above_ground':cmin-e0,
        'retained_weight_gap_upper':gap_upper,
        'response_norm_ground_lower':response_norm_floor,
        'shifted_target_energy':energy,'shifted_gap':delta,
        'shifted_spectral_upper':upper,'condition_number':condition,
        'coupling_operator_norm':g,'chebyshev_degree_for_0p001t':chebyshev_degree,
        'diagonal_relative_residual_norm':rho,
        'first_sampled_jacobi_iteration_below_margin':first,
        'jacobi_residual_history':history,
        'matrix_nonzeros':int(h.nnz),'eliminated_nonzeros':int(c.nnz),
        'seconds':time.monotonic()-start,
        'status':'enumerated_float64_diagnostic_not_certified'
    }


def case(rungs, coupling, max_steps=32, margin=0.001):
    start=time.monotonic()
    labels,rows=rational_rows(rungs,F(coupling))
    h=numerical(rows)
    e0,psi=smallest(h,vector=True)
    pidx=np.array([i for i,v in enumerate(labels) if doublons(v,2*rungs)==0])
    bare=analyze(h,pidx,psi,e0,margin,max_steps)
    rotated=h; psi_rotated=psi.copy();gate=local_rotation()
    for rung in range(rungs):
        G=lifted_gate(labels,rung,gate)
        rotated=(G.T@rotated@G).tocsr();rotated.eliminate_zeros()
        psi_rotated=G.T@psi_rotated
    dressed=analyze(rotated,pidx,psi_rotated,e0,margin,max_steps)
    return {
        'rungs':rungs,'sites':2*rungs,'U_over_t':8,
        'horizontal_hopping_over_t':str(F(coupling)),
        'vertical_hopping_over_t':'1','boundary':'open','filling':'half',
        'N_up':rungs,'N_down':rungs,'global_determinants_constructed':len(labels),
        'ground_energy_over_t':e0,
        'rational_local_dressing_half_angle':str(RATIONAL_HALF_ANGLE),
        'dressing_uses_global_ground_state':False,
        'state_used_only_for_charged_diagnostic':True,
        'bare':bare,'rational_dimer_dressed':dressed,
        'seconds':time.monotonic()-start,
    }


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--rungs',type=int,nargs='+',required=True)
    p.add_argument('--couplings',nargs='+',default=['0','1/4','1/2','1'])
    p.add_argument('--steps',type=int,default=32)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if any(L not in (1,2,3,4) for L in a.rungs) or not 0 <= a.steps <= 64:
        raise ValueError('Bounded finite diagnostic only')
    if a.out.exists():
        raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    for L in a.rungs:
        for lam in a.couplings:
            result=case(L,lam,a.steps)
            path=a.out/('ladder_%s_lambda_%s.json'%(L,str(F(lam)).replace('/','_')))
            path.write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({'case':path.name,'seconds':result['seconds'],
                'E0':result['ground_energy_over_t'],
                'bare_weight':result['bare']['retained_ground_weight'],
                'dressed_weight':result['rational_dimer_dressed']['retained_ground_weight'],
                'bare_gap':result['bare']['eliminated_gap_above_ground'],
                'dressed_gap':result['rational_dimer_dressed']['eliminated_gap_above_ground']}),flush=True)
