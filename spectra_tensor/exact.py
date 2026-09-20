"""Independent exact acceptance for finite Hubbard scalar responses.

No numerical package is imported. The optional compiled backend performs only
Gaussian-integer contractions; the pure-Python backend implements the same
arithmetic for cross-checks. Neither evaluates a full determinant basis.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import lcm
from pathlib import Path
from time import perf_counter
import ctypes as ct
import json, os, resource, subprocess, sys, tempfile, shlex

PHYS=((0,0),(1,0),(0,1),(1,1))
def plus(a,b): return (a[0]+b[0],a[1]+b[1])
def rat(x):
    if type(x) not in (str,int,F): raise ValueError('literal rational required')
    return F(x)
def integer(x):
    if type(x) is not int: raise ValueError('literal integer required')
    return x

def validate_model(spec):
    if spec.get('kind')!='spin_independent_hubbard_v1': raise ValueError('model kind')
    n,p=spec.get('sites'),spec.get('particles')
    if type(n) is not int or not 2<=n<=200 or n%2 or type(p) is not int or p!=n:
        raise ValueError('even half-filled finite model required')
    if spec.get('mode_order')!='up_then_down' or spec.get('units')!='t': raise ValueError('convention')
    u=[rat(x) for x in spec['U']]
    if len(u)!=n or any(x<0 for x in u): raise ValueError('repulsive onsite terms')
    edges=[];seen=set()
    for a,b,v in spec['edges']:
        integer(a);integer(b)
        if not 0<=a<b<n or (a,b) in seen: raise ValueError('edge convention')
        seen.add((a,b));edges.append((a,b,rat(v)))
    return n,u,edges

def model_hash(spec):
    validate_model(spec)
    return sha256(json.dumps({k:v for k,v in spec.items() if k!='tag'},sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()).hexdigest()

def source(spec): return [1 if i%2==0 else 2 for i in range(validate_model(spec)[0])]

def validate_program(program,n,max_bond=256):
    if program.get('kind')!='charge_mps_gaussian_integer_v1': raise ValueError('program kind')
    qs,cores=program['charges'],program['cores'];scale=integer(program['local_denominator'])
    if not 1<=scale<=1<<40: raise ValueError('denominator envelope')
    if len(qs)!=n+1 or len(cores)!=n or qs[0]!=[[0,0]] or qs[-1]!=[[n//2,n//2]]: raise ValueError('boundary sector')
    for cut,charges in enumerate(qs):
        if not 1<=len(charges)<=max_bond: raise ValueError('bond envelope')
        for q in charges:
            if len(q)!=2 or any(type(v) is not int or not 0<=v<=min(cut,n//2) for v in q): raise ValueError('bond charge')
    for i,core in enumerate(cores):
        seen=set()
        for entry in core:
            if len(entry)!=5 or any(type(v) is not int for v in entry): raise ValueError('integer core entry')
            l,p,r,re,im=entry
            if not 0<=l<len(qs[i]) or not 0<=p<4 or not 0<=r<len(qs[i+1]): raise ValueError('core index')
            if (l,p,r) in seen: raise ValueError('duplicate entry')
            seen.add((l,p,r))
            if max(abs(re),abs(im))>=1<<48 or not (re or im): raise ValueError('coefficient envelope')
            if plus(qs[i][l],PHYS[p])!=tuple(qs[i+1][r]): raise ValueError('particle sector leakage')
    return program

def integer_hamiltonian(spec):
    """Exact interleaved Jordan-Wigner automaton for den*H.

    Every accepting path uses exactly one rational coefficient. Consequently
    den is global, not raised to the number of sites.
    """
    n,u,edges=validate_model(spec);den=lcm(*(v.denominator for v in u),*(v.denominator for _,_,v in edges))
    unit=[(j,j,1) for j in range(4)];parity=[(j,j,(-1)**sum(PHYS[j])) for j in range(4)]
    ann=[[(q^(1<<s),q,(-1)**((q&((1<<s)-1)).bit_count())) for q in range(4) if (q>>s)&1] for s in (0,1)]
    pending=[(ei,a,b,s,sgn,int(v*den)) for ei,(a,b,v) in enumerate(edges) if v for s in (0,1) for sgn in (1,-1)]
    labels=[];charges=[]
    for cut in range(n+1):
        ls=['start'] if cut==0 else ['end'] if cut==n else ['start','end']+[(ei,s,sgn) for ei,a,b,s,sgn,v in pending if a<cut<=b]
        labels.append(ls);charges.append([(0,0) if isinstance(k,str) else (k[2],0) if k[1]==0 else (0,k[2]) for k in ls])
    cores=[]
    for i in range(n):
        left={v:k for k,v in enumerate(labels[i])};right={v:k for k,v in enumerate(labels[i+1])};ops=[]
        def put(a,b,op):
            if a in left and b in right: ops.extend((left[a],right[b],p,q,v) for p,q,v in op if v)
        put('start','start',unit);put('end','end',unit);put('start','end',[(3,3,int(u[i]*den))])
        for ei,a,b,s,sgn,v in pending:
            key=(ei,s,sgn)
            if i==a:
                mat=[(q,p,k*((-1)**sum(PHYS[p]))) for p,q,k in ann[s]] if sgn==1 else [(p,q,((-1)**sum(PHYS[p]))*k) for p,q,k in ann[s]]
                put('start',key,[(p,q,-v*k) for p,q,k in mat])
            elif a<i<b: put(key,key,parity)
            elif i==b: put(key,'end',ann[s] if sgn==1 else [(q,p,k) for p,q,k in ann[s]])
        cores.append(ops)
    return dict(cores=cores,charges=charges,den=den)

class ActionCores:
    """Reiterable local action: the expanded chain is never stored."""
    def __init__(self,mpo,program):
        self.mpo,self.program=mpo,program;self.counts=[None]*len(program['cores'])
    def __len__(self): return len(self.counts)
    def __getitem__(self,i):
        if not 0<=i<len(self): raise IndexError(i)
        qs=self.program['charges'];grouped=defaultdict(list)
        for l,p,r,re,im in self.program['cores'][i]: grouped[p].append((l,r,re,im))
        new={};dl,dr=len(qs[i]),len(qs[i+1])
        for a,b,p,q,v in self.mpo['cores'][i]:
            for l,r,re,im in grouped[q]:
                key=(a*dl+l,p,b*dr+r);u,w=new.get(key,(0,0));new[key]=(u+v*re,w+v*im)
        entries=[[*key,re,im] for key,(re,im) in sorted(new.items()) if re or im]
        if any(max(abs(e[3]),abs(e[4]))>=1<<60 for e in entries): raise ValueError('action coefficient envelope')
        self.counts[i]=len(entries);return entries
    def __iter__(self):
        for i in range(len(self)): yield self[i]
    @property
    def entry_count(self):
        for i,count in enumerate(self.counts):
            if count is None: self[i]
        return sum(self.counts)

def apply_integer(mpo,program,lazy=True):
    qs=[[plus(q,d) for d in delta for q in old] for old,delta in zip(program['charges'],mpo['charges'])]
    cores=ActionCores(mpo,program)
    return dict(charges=qs,cores=cores if lazy else list(cores),local_denominator=program['local_denominator'])

def iter_blocks(program):
    def group(charges):
        g=defaultdict(list)
        for i,q in enumerate(charges): g[tuple(q)].append(i)
        return g,{index:j for indices in g.values() for j,index in enumerate(indices)}
    left,lp=group(program['charges'][0])
    for i,entries in enumerate(program['cores']):
        right,rp=group(program['charges'][i+1]);current={}
        for l,p,r,re,im in entries:
            q=tuple(program['charges'][i][l]);key=(q,p);qr=plus(q,PHYS[p])
            if key not in current:
                nr,nc=len(left[q]),len(right[qr]);current[key]=(nr,nc,[0]*(nr*nc),[0]*(nr*nc))
            nr,nc,a,b=current[key];k=lp[l]*nc+rp[r];a[k]=re;b[k]=im
        yield current
        left,lp=right,rp

class Native:
    def __init__(self):
        src=Path(__file__).with_name('exact_contract.c');self.source_hash=sha256(src.read_bytes()).hexdigest()
        # Private per-user cache. Existing source-hash names identify the build.
        root=Path(tempfile.gettempdir())/f'spectra_exact_{os.getuid()}';root.mkdir(mode=0o700,exist_ok=True)
        target=root/f'contract_{self.source_hash}.so'
        if not target.exists():
            tmp=root/f'build_{os.getpid()}.so'
            cmd=shlex.split(os.environ.get('CC','cc'))+['-O3','-std=c11','-shared','-fPIC']
            prefix=os.environ.get('SPECTRA_GMP_PREFIX')
            if prefix:cmd.extend(['-I'+str(Path(prefix)/'include'),'-L'+str(Path(prefix)/'lib')])
            cmd+=shlex.split(os.environ.get('CPPFLAGS',''))+[str(src),'-o',str(tmp)]+shlex.split(os.environ.get('LDFLAGS',''))+['-lgmp']
            r=subprocess.run(cmd,capture_output=True,text=True)
            if r.returncode: raise RuntimeError('GMP compilation failed: '+r.stderr)
            tmp.replace(target)
        lib=self.lib=ct.CDLL(str(target));P=ct.c_void_p;A=ct.POINTER(ct.c_int64);S=ct.c_size_t
        lib.matrix_one.argtypes=[];lib.matrix_one.restype=P
        lib.matrix_transfer.argtypes=[P,A,A,S,A,A,S,ct.c_int];lib.matrix_transfer.restype=P
        lib.matrix_add.argtypes=[P,P];lib.matrix_add.restype=ct.c_int
        lib.matrix_free.argtypes=[P];lib.matrix_free.restype=None
        lib.matrix_scalar.argtypes=[P,ct.c_int];lib.matrix_scalar.restype=P
        lib.string_free.argtypes=[P];lib.string_free.restype=None
    def scalar(self,p):
        values=[]
        for imaginary in (0,1):
            text=self.lib.matrix_scalar(p,imaginary)
            if not text: raise RuntimeError('non-scalar boundary')
            try: values.append(int(ct.string_at(text)))
            finally: self.lib.string_free(text)
        return tuple(values)

def exact_overlap(a,b,native=None,conjugate=True):
    pairs=((v,v) for v in iter_blocks(a)) if a is b else zip(iter_blocks(a),iter_blocks(b))
    if native is None: return python_overlap(pairs,conjugate)
    lib=native.lib;env={(0,0):lib.matrix_one()};ops=0
    if not env[(0,0)]: raise MemoryError('exact boundary allocation')
    try:
        for ax,bx in pairs:
            old=env;env={}
            try:
                for q,e in old.items():
                    for p,ph in enumerate(PHYS):
                        aa,bb=ax.get((q,p)),bx.get((q,p))
                        if aa is None or bb is None: continue
                        m,ra,ar,ai=aa;n,rb,br,bi=bb
                        arrays=[(ct.c_int64*len(v))(*v) for v in (ar,ai,br,bi)]
                        term=lib.matrix_transfer(e,arrays[0],arrays[1],ra,arrays[2],arrays[3],rb,int(conjugate))
                        if not term: raise MemoryError('exact transfer allocation')
                        qnew=plus(q,ph)
                        if qnew in env:
                            ok=lib.matrix_add(env[qnew],term);lib.matrix_free(term)
                            if not ok: raise AssertionError('exact environment shape')
                        else: env[qnew]=term
                        ops+=1
            finally:
                for ptr in old.values(): lib.matrix_free(ptr)
        if not env: return (0,0),ops
        if len(env)!=1: raise AssertionError('non-scalar boundary')
        return native.scalar(next(iter(env.values()))),ops
    finally:
        for ptr in env.values(): lib.matrix_free(ptr)

def python_overlap(pairs,conjugate):
    env={(0,0):(1,1,[(1,0)])};ops=0
    for ax,bx in pairs:
        out={}
        for q,(m,n,e) in env.items():
            for p,ph in enumerate(PHYS):
                aa,bb=ax.get((q,p)),bx.get((q,p))
                if aa is None or bb is None: continue
                _,ra,ar,ai=aa;_,rb,br,bi=bb;temp=[(0,0)]*(m*rb)
                for i in range(m):
                    for j in range(rb):
                        re=im=0
                        for k in range(n):
                            er,ei=e[i*n+k];vr,vi=br[k*rb+j],bi[k*rb+j];re+=er*vr-ei*vi;im+=er*vi+ei*vr
                        temp[i*rb+j]=(re,im)
                qnew=plus(q,ph)
                if qnew not in out: out[qnew]=(ra,rb,[(0,0)]*(ra*rb))
                dest=out[qnew][2]
                for i in range(ra):
                    for j in range(rb):
                        re,im=dest[i*rb+j]
                        for k in range(m):
                            tr,ti=temp[k*rb+j];vr=ar[k*ra+i];vi=ai[k*ra+i]*(-1 if conjugate else 1);re+=vr*tr-vi*ti;im+=vr*ti+vi*tr
                        dest[i*rb+j]=(re,im)
                ops+=1
        env=out
    return (next(iter(env.values()))[2][0] if env else (0,0)),ops

def source_amplitude(program,physical):
    v={0:(1,0)}
    for core,p in zip(program['cores'],physical):
        new={}
        for l,q,r,re,im in core:
            if q!=p or l not in v: continue
            a,b=v[l];u,w=new.get(r,(0,0));new[r]=(u+a*re-b*im,w+a*im+b*re)
        v=new
    return v.get(0,(0,0))

def verify(candidate,request,backend='gmp'):
    start=perf_counter()
    if request.get('kind')!='tensor_hubbard_response_request_v1' or candidate.get('kind')!='tensor_hubbard_response_v1': raise ValueError('request/candidate kind')
    spec=request['model'];n,_,_=validate_model(spec);identity=model_hash(spec)
    if candidate.get('model_hash')!=identity: raise ValueError('requested model mismatch')
    physical=source(spec)
    if request.get('source')!=physical: raise ValueError('requested source')
    expected=request['frequencies'];queries=candidate['queries']
    if not expected or len(expected)!=len(queries): raise ValueError('frequency count')
    if backend not in ('python','gmp'): raise ValueError('arithmetic backend')
    native=Native() if backend=='gmp' else None;mpo=integer_hamiltonian(spec);hd=mpo['den'];results=[];transfers=0
    for frequency,query in zip(expected,queries):
        w,eta=rat(frequency['omega']),rat(frequency['eta'])
        if eta<=0 or any(rat(query[k])!=rat(frequency[k]) for k in ('omega','eta')): raise ValueError('requested frequency')
        x=validate_program(query['program'],n);hx=apply_integer(mpo,x);values={}
        for name,a,b,conj in [('s',x,x,True),('k',x,hx,True),('l',hx,hx,True),('xx',x,x,False),('xh',x,hx,False)]:
            values[name],ops=exact_overlap(a,b,native,conj);transfers+=ops
        if any(values[k][1] for k in ('s','k','l')): raise AssertionError('nonreal Hermitian moment')
        if values['s'][0]<0 or values['l'][0]<0: raise AssertionError('negative norm')
        D=x['local_denominator']**n;D2=D*D;bxr,bxi=source_amplitude(x,physical);bhr,bhi=source_amplitude(hx,physical)
        s=F(values['s'][0],D2);k=F(values['k'][0],D2*hd);ell=F(values['l'][0],D2*hd*hd)
        r2=1-2*(w*F(bxr,D)-eta*F(bxi,D)-F(bhr,D*hd))+(w*w+eta*eta)*s-2*w*k+ell
        if r2<0: raise AssertionError('negative exact residual norm')
        xxr,xxi=values['xx'];xhr,xhi=values['xh']
        re=2*F(bxr,D)-w*F(xxr,D2)+eta*F(xxi,D2)+F(xhr,D2*hd)
        im=2*F(bxi,D)-w*F(xxi,D2)-eta*F(xxr,D2)+F(xhi,D2*hd)
        radius=r2/eta;lo=max(F(0),-im-radius);hi=min(1/eta,-im+radius)
        if lo>hi: raise AssertionError('empty spectral interval')
        results.append(dict(omega=str(w),eta=str(eta),center_real=str(re),center_imag=str(im),radius=str(radius),residual_norm_squared=str(r2),absorption_lower=str(lo),absorption_upper=str(hi),max_bond=max(map(len,x['charges'])),integer_entries=sum(map(len,x['cores'])),action_entries=hx['cores'].entry_count))
    return dict(status='accepted_exact_tensor_response',model_hash=identity,queries=results,global_sector_enumerated=False,determinants_enumerated=0,arithmetic='Gaussian integers and rationals',backend=backend,native_source_sha256=native.source_hash if native else None,verifier_source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),block_transfers=transfers,seconds=perf_counter()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),numerical_packages_loaded=[k for k in ('numpy','scipy','torch') if k in sys.modules],scope='Declared finite Hubbard model and source, positive-broadening scalar resolvent; no ground-state or continuum claim.')

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('candidate',type=Path);p.add_argument('--request',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--backend',choices=['gmp','python'],default='gmp');a=p.parse_args()
    if a.out.exists(): raise FileExistsError(a.out)
    craw=a.candidate.read_bytes();rraw=a.request.read_bytes();r=verify(json.loads(craw),json.loads(rraw),a.backend)
    r.update(candidate_sha256=sha256(craw).hexdigest(),request_sha256=sha256(rraw).hexdigest());a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x') as f: json.dump(r,f,indent=2)
    print(json.dumps({k:v for k,v in r.items() if k!='queries'},indent=2))
if __name__=='__main__': main()
