"""Numerical response proposal and residual screening; acceptance is separate."""
from time import perf_counter
import numpy as np
from . import tensor as t,physics as p,exact,io

def screen(spec,x,z,bits=30):
    if not np.isfinite(z) or z.imag<=0:raise ValueError('finite positive broadening')
    program=io.export_program(x,bits);h=exact.integer_hamiltonian(spec);hp=exact.apply_integer(h,program);scale=program['local_denominator']
    def prepare(prog):
        for core in exact.iter_blocks(prog):
            yield {key:(np.array(a,float)+1j*np.array(b,float)).reshape(nr,nc)/scale for key,(nr,nc,a,b) in core.items()}
    def dot(a,b,conjugate=True):
        env={(0,0):np.ones((1,1),complex)}
        pairs=((core,core) for core in prepare(a)) if a is b else zip(prepare(a),prepare(b))
        for ax,bx in pairs:
            out={}
            for q,e in env.items():
                for pp,ph in enumerate(t.PHYS):
                    A,B=ax.get((q,pp)),bx.get((q,pp))
                    if A is None or B is None:continue
                    term=(A.conj().T if conjugate else A.T)@e@B;qn=t.plus(q,ph)
                    if qn in out:out[qn]+=term
                    else:out[qn]=term
            env=out
        return next(iter(env.values()))[0,0] if env else 0.
    D=scale**len(program['cores']);hd=h['den'];ar,ai=exact.source_amplitude(program,p.source(spec));br,bi=exact.source_amplitude(hp,p.source(spec))
    bx=ar/D+1j*(ai/D);bh=br/(D*hd)+1j*(bi/(D*hd))
    s=dot(program,program).real;k=dot(program,hp).real/hd;l=dot(hp,hp).real/hd**2
    r2=float(1-2*(z*bx-bh).real+abs(z)**2*s-2*z.real*k+l)
    if not np.isfinite(r2) or r2<-1e-8:raise ArithmeticError('unreliable numerical residual')
    center=2*bx-z*dot(program,program,False)+dot(program,hp,False)/hd
    return dict(radius=max(0.,r2)/z.imag,residual_squared=r2,center_real=float(center.real),center_imag=float(center.imag),bond=x.bond,entries=x.entries)

def cocg(spec,z,max_bond=48,maxiter=12,rtol=.02,cutoff=1e-10,assess_every=12,initial=None,verbose=True):
    if maxiter<1 or assess_every<1 or rtol<=0:raise ValueError('iteration parameters')
    start=perf_counter();mpo=p.compile_hubbard(spec);b=t.product(p.source(spec));x=b.scaled(1/z) if initial is None else initial.copy();history=[];best=None;discard=0.
    def rounded(a):
        nonlocal discard
        try:a,d=t.compress(a,max_bond,cutoff)
        except ArithmeticError:
            if t.norm(a)>1e-30:raise
            return b.scaled(0.)
        discard+=d;return a
    result=screen(spec,x,z)
    if result['radius']<=rtol**2/z.imag:return x,dict(history=[dict(step=0,**result)],seconds=perf_counter()-start,enumerated_states=0)
    hx,_=t.apply_compressed(mpo,x,max_bond,cutoff);r=rounded(t.add(t.add(b,hx),x,beta=-z));v=r.copy();rho=t.overlap(r,r,False)
    for step in range(1,maxiter+1):
        hv,_=t.apply_compressed(mpo,v,max_bond,cutoff);av=rounded(t.add(v,hv,alpha=z,beta=-1))
        den=t.overlap(v,av,False)
        if abs(den)<1e-28 or abs(rho)<1e-28:break
        alpha=rho/den;x=rounded(t.add(x,v,beta=alpha));r=rounded(t.add(r,av,beta=-alpha));residual=t.norm(r)
        if step%assess_every==0 or residual<rtol or step==maxiter:
            result=screen(spec,x,z);row=dict(step=step,recurrence_residual=residual,radius=result['radius'],bond=x.bond,seconds=perf_counter()-start);history.append(row)
            if verbose:print(row,flush=True)
            if best is None or row['radius']<best[0]:best=(row['radius'],x.copy())
            if row['radius']<=rtol**2/z.imag:break
            if residual<rtol or step%20==0:
                hx,_=t.apply_compressed(mpo,x,max_bond,cutoff);r=rounded(t.add(t.add(b,hx),x,beta=-z));v=r.copy();rho=t.overlap(r,r,False);row['restart']=True;continue
        rn=t.overlap(r,r,False);v=rounded(t.add(r,v,beta=rn/rho));rho=rn
    if best is None:best=(screen(spec,x,z)['radius'],x.copy())
    return best[1],dict(history=history,seconds=perf_counter()-start,cumulative_discard_diagnostic=discard,enumerated_states=0)
