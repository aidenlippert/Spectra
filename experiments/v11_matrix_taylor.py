"""Conventional exact matrix implementation of fraction-free Taylor action.

No Liouville matrix, floating coefficients, learned method or new certificate.
Fixed-width kernels run only after integer forward overflow bounds pass;
otherwise NumPy object arrays carry exact Python integers. Output is replayed
by the original V7 residual checker by the caller.
"""
from fractions import Fraction as F
from math import factorial,lcm
import numpy as np
from .v7_certificate import Generator,Piece,clean,rational,norm_witness,BudgetExceeded

INT_MAX=(1<<63)-1


def bounded_integer(z):
    z=int(z)
    if abs(z).bit_length()>8192:raise BudgetExceeded('matrix integer bit budget')
    return z


def masks(p):
    x=z=0
    for c in p:
        x=(x<<1)|int(c in 'XY');z=(z<<1)|int(c in 'YZ')
    return x,z,(x&z).bit_count()


def label(x,z,n):
    return ''.join('IXZY'[(x>>i&1)+2*(z>>i&1)] for i in reversed(range(n)))


def rotate(r,im,phase):
    phase%=4
    if phase==0:return r,im
    if phase==1:return -im,r
    if phase==2:return -r,-im
    return im,-r


def magnitude(pair):
    # Avoid abs(int64_min), which itself overflows. Stored arrays never contain
    # that value, but taking min/max then converting to Python also checks safely.
    return max(abs(int(v)) for a in pair for v in (a.min(),a.max()))


class MatrixAction:
    def __init__(self,h,gamma,n,cap=512):
        if type(n)is not int or not 1<=n<=6:raise BudgetExceeded('matrix n<=6 budget')
        g=Generator(h,gamma,n,cap)
        self.n,self.d,self.cap=n,1<<n,cap
        self.h,self.gamma=g.h,g.gamma
        self.q=1
        for c in [*self.h.values(),self.gamma/4]:
            self.q=bounded_integer(lcm(self.q,c.denominator))
        self.hi={p:bounded_integer(c*self.q) for p,c in self.h.items() if p!='I'*n}
        self.gi=bounded_integer(self.gamma*self.q/4)
        self.index=np.arange(self.d,dtype=np.int64)
        self.cost=dict(matrix_entries=2*self.d*self.d,initial_updates=0,
                       action_calls=0,action_entry_updates=0,permutation_builds=0,
                       transforms=0,walsh_add_subtracts=0,exported_integers=0,
                       int64_kernels=0,object_kernels=0,array_cast_entries=0,
                       range_scan_entries=0,max_integer_bits=0,max_forward_bits=0,
                       peak_pauli_terms=0)
        self.perms={}
        for p in [*self.hi,*[('I'*i+a+'I'*(n-i-1)) for i in range(n) for a in 'XYZ']]:
            if p in self.perms:continue
            x,z,y=masks(p)
            signs=np.array([1 if not (int(s)&z).bit_count()%2 else -1 for s in self.index],dtype=np.int64)
            perm=self.index^x
            self.perms[p]=(perm,signs,signs[perm],y)
            self.cost['permutation_builds']+=1

    def checked_arrays(self,pair,bound):
        bound=bounded_integer(bound)
        self.cost['max_forward_bits']=max(self.cost['max_forward_bits'],bound.bit_length())
        dtype=np.int64 if bound<=INT_MAX else object
        key='int64_kernels' if dtype is np.int64 else 'object_kernels'
        self.cost[key]+=1
        result=[]
        for a in pair:
            if a.shape!=(self.d,self.d) or a.dtype not in (np.dtype('int64'),np.dtype('O')):raise ValueError('integer matrix shape/dtype')
            if a.dtype!=np.dtype(dtype):self.cost['array_cast_entries']+=a.size
            result.append(a.astype(dtype,copy=False))
        return tuple(result)

    def initial(self,op):
        op=clean(op,self.n,self.cap);q0=1
        for c in op.values():q0=bounded_integer(lcm(q0,c.denominator))
        ints={p:bounded_integer(c*q0) for p,c in op.items()}
        bound=bounded_integer(sum(map(abs,ints.values())))
        dtype=np.int64 if bound<=INT_MAX else object
        r=np.zeros((self.d,self.d),dtype=dtype);im=r.copy()
        for p,c in ints.items():
            x,z,y=masks(p)
            s=np.array([1 if not (int(j)&z).bit_count()%2 else -1 for j in self.index],dtype=dtype)
            if y%2:
                im[self.index^x,self.index]+=c*s*(1 if y%4==1 else -1)
            else:r[self.index^x,self.index]+=c*s*(1 if y%4==0 else -1)
            self.cost['initial_updates']+=self.d
        self.cost['max_integer_bits']=bound.bit_length()
        return (r,im),q0

    def apply(self,pair):
        m=magnitude(pair);self.cost['range_scan_entries']+=2*self.d*self.d
        bound=bounded_integer((2*sum(map(abs,self.hi.values()))+6*self.n*abs(self.gi))*m)
        r,im=self.checked_arrays(pair,max(m,bound))
        out_r=np.zeros_like(r);out_i=np.zeros_like(im)
        if not m:
            self.cost['action_calls']+=1
            return out_r,out_i
        for p,c in self.hi.items():
            perm,sign,signperm,y=self.perms[p]
            dr=r[perm,:]*signperm[:,None]-r[:,perm]*sign[None,:]
            di=im[perm,:]*signperm[:,None]-im[:,perm]*sign[None,:]
            dr,di=rotate(dr,di,y+1)
            out_r+=c*dr;out_i+=c*di
            self.cost['action_entry_updates']+=2*self.d*self.d
        if self.gi:
            for i in range(self.n):
                for a in 'XYZ':
                    p='I'*i+a+'I'*(self.n-i-1)
                    perm,sign,signperm,y=self.perms[p]
                    factor=(-1 if y%2 else 1)*signperm[:,None]*sign[None,:]
                    dr=r[perm[:,None],perm[None,:]]*factor-r
                    di=im[perm[:,None],perm[None,:]]*factor-im
                    out_r+=self.gi*dr;out_i+=self.gi*di
                    self.cost['action_entry_updates']+=2*self.d*self.d
        self.cost['action_calls']+=1
        actual=magnitude((out_r,out_i));self.cost['range_scan_entries']+=2*self.d*self.d
        self.cost['max_integer_bits']=max(self.cost['max_integer_bits'],actual.bit_length())
        return out_r,out_i

    def to_pauli(self,pair):
        m=magnitude(pair);self.cost['range_scan_entries']+=2*self.d*self.d
        r,im=self.checked_arrays(pair,bounded_integer(self.d*m))
        # Each row x contains A[s,s xor x], s varying over columns.
        s=self.index[None,:];x=self.index[:,None]
        arrays=[a[s,s^x].copy() for a in (r,im)]
        half=1
        while half<self.d:
            for a in arrays:
                blocks=a.reshape(self.d,-1,2,half)
                left=blocks[:,:,0,:].copy();right=blocks[:,:,1,:].copy()
                blocks[:,:,0,:]=left+right;blocks[:,:,1,:]=left-right
            self.cost['walsh_add_subtracts']+=2*self.d*self.d
            half*=2
        out={}
        for x in range(self.d):
            for z in range(self.d):
                re,imag=rotate(int(arrays[0][x,z]),int(arrays[1][x,z]),(x&z).bit_count())
                if imag:raise ValueError('non-Hermitian exact matrix')
                if re%self.d:raise ValueError('nonintegral Pauli-coordinate matrix')
                if re:out[label(x,z,self.n)]=bounded_integer(re//self.d)
        if len(out)>self.cap:raise BudgetExceeded('matrix Pauli support budget')
        self.cost['transforms']+=1;self.cost['exported_integers']+=len(out)
        self.cost['peak_pauli_terms']=max(self.cost['peak_pauli_terms'],len(out))
        return out


def matrix_taylor(gen,initial,duration,tolerance,max_order=24):
    if type(max_order)is not int or not 0<=max_order<=24:raise ValueError('Taylor order cap')
    duration,tolerance=rational(duration),rational(tolerance)
    if duration<=0 or tolerance<0:raise ValueError('duration/tolerance')
    initial=clean(initial,gen.n,gen.max_terms)
    engine=MatrixAction(gen.h,gen.gamma,gen.n,gen.max_terms)
    b,q0=engine.initial(initial);raw=engine.to_pauli(b)
    coeffs=[];probes=[];pairs=sorts=0
    try:
        for k in range(max_order+1):
            denom=bounded_integer(q0*engine.q**k*factorial(k))
            cs={p:F(c,denom) for p,c in raw.items()};coeffs.append(cs)
            b=engine.apply(b);raw=engine.to_pauli(b)
            derivative_denom=bounded_integer(denom*engine.q)
            derivative={p:F(c,derivative_denom) for p,c in raw.items()}
            for grouping in ('l1','firstfit','weighted'):
                groups,cost=norm_witness(derivative,grouping)
                pairs+=cost['group_comparisons'];sorts+=cost['sorted_terms']
                bound=duration**(k+1)/F(k+1)*sum((F(g['upper']) for g in groups),F(0))
                probes.append(dict(order=k,grouping=grouping,bound=str(bound)))
                if bound<=tolerance:
                    witnesses={'jump:0':[],**{f'residual:0:{j}':[] for j in range(k+1)}}
                    witnesses[f'residual:0:{k}']=groups if derivative else []
                    if not derivative:witnesses={'jump:0':[],'residual:0:0':[]}
                    w=dict(schema='v7-residual-1',integration_basis='power',witnesses=witnesses,
                           claimed_bound=str(bound),construction_cost=dict(group_comparisons=pairs,sorted_terms=sorts))
                    return Piece(duration,tuple(coeffs)),w,dict(probes=probes,matrix_work=dict(engine.cost),denominator=engine.q,initial_denominator=q0)
        raise BudgetExceeded('matrix Taylor order without certificate')
    except ValueError as exc:
        exc.construction_cost=dict(probes=probes,matrix_work=dict(engine.cost),denominator=engine.q,initial_denominator=q0)
        raise
