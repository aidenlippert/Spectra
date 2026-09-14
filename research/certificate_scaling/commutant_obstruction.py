"""Exact inertia certificate excluding near-density quartic interactions.

See approximate_commutant_obstruction.md. Real rational input; optionally all
complex orbital unitaries and one-body number-ideal multipliers. The norm is
the two-particle quartic coefficient Frobenius norm, not a claimed
lower bound on fixed-N operator norm or ground-state energy error.
"""
from fractions import Fraction as F
from math import lcm,isqrt
from pathlib import Path
import argparse,hashlib,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import decode,number_shift,mono,product,add,scale,hermitian
from research.certificate_scaling.hidden_density_basis import commutator_map,generators


def number_ideal_core(h,m):
    """Orthogonally remove quartic parts of Nhat*Q(A), for all one-body A.

    For real Hermitian input the imaginary-antisymmetric lift is orthogonal.
    This projection is covariant under complex orbital unitaries.
    """
    if m<=2:raise ValueError('Number-ideal projection requires M>2')
    if any(not isinstance(c,(int,F)) for c in h.values()) or not hermitian(h):
        raise ValueError('Projection implementation requires real rational Hermitian input')
    v={w:c for w,c in h.items() if len(w)==4};lifts=[];moments=[]
    for i,j in generators(m):
        q=mono(((1,i),(0,j))) if i==j else add(mono(((1,i),(0,j))),mono(((1,j),(0,i))))
        lift={w:c for w,c in product(number_shift(m,0),q).items() if len(w)==4}
        lifts.append(lift);moments.append(sum(c*v.get(w,F(0)) for w,c in lift.items()))
    diagonal_sum=sum(t for t,(i,j) in zip(moments,generators(m)) if i==j)
    coefficients=[(t-diagonal_sum/F(2*(m-1)))/(m-2) if i==j else t/F(2*(m-2))
                  for t,(i,j) in zip(moments,generators(m))]
    core=add(v,*(scale(lift,-c) for lift,c in zip(lifts,coefficients)))
    if any(sum(c*core.get(w,F(0)) for w,c in lift.items()) for lift in lifts):
        raise ArithmeticError('Projection not orthogonal to one-body lift')
    # These are the remaining Hermitian generators after multiplication by i.
    # Their projection coefficients vanish for real Hermitian V, but check that
    # fact exactly rather than silently relying on the real-input convention.
    for i in range(m):
        for j in range(i+1,m):
            q=add(mono(((1,i),(0,j))),mono(((1,j),(0,i)),-1))
            lift={w:c for w,c in product(number_shift(m,0),q).items() if len(w)==4}
            if sum(c*core.get(w,F(0)) for w,c in lift.items()):
                raise ArithmeticError('Core not orthogonal to imaginary one-body lift')
    return core


def integer_gram(h,m,antisymmetric=False):
    columns,words=commutator_map(h,m,antisymmetric)
    denominator=lcm(*(c.denominator for col in columns for c in col.values()))
    rows={w:[] for w in words}
    for j,col in enumerate(columns):
        for w,c in col.items():rows[w].append((j,int(c*denominator)))
    n=len(columns);g=[[0]*n for _ in range(n)]
    for entries in rows.values():
        for k,(i,x) in enumerate(entries):
            for j,y in entries[k:]:g[i][j]+=x*y
    for i in range(n):
        for j in range(i):g[i][j]=g[j][i]
    return g,denominator,len(words),sum(map(len,columns))


def inertia_integer(matrix):
    """Fraction-free symmetric elimination, exact congruence signs.

    Refuses a zero leading pivot; it never guesses inertia at a degeneracy.
    With nonzero pivots, Bareiss pivots are leading principal determinants;
    LDL pivots are successive ratios, whose signs give inertia.
    """
    n=len(matrix)
    if any(len(row)!=n for row in matrix) or any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(i)):
        raise ValueError('Symmetric square integer matrix required')
    if any(type(x) is not int for row in matrix for x in row):raise ValueError('Integer matrix required')
    a=[row[:] for row in matrix];previous=1;positive=negative=0;pivots=[];peak=0
    for k in range(n):
        pivot=a[k][k]
        if not pivot:raise ValueError('Zero leading pivot: choose another rational shift')
        positive+=pivot*previous>0;negative+=pivot*previous<0;pivots.append(pivot)
        for i in range(k+1,n):
            for j in range(i,n):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                value,remainder=divmod(numerator,previous)
                if remainder:raise ArithmeticError('Non-exact Bareiss division')
                a[i][j]=a[j][i]=value
                peak=max(peak,abs(value).bit_length())
        previous=pivot
    digest=hashlib.sha256(','.join(hex(x) for x in pivots).encode()).hexdigest()
    return {'positive':positive,'negative':negative,'zero':0,'leading_pivot_sha256':digest,'peak_integer_bits':peak}


def obstruction(h,m,threshold=F(1,10000),complex_unitary=False,quotient_number_ideal=False):
    start=time.monotonic()
    if threshold<=0:raise ValueError('Positive threshold required')
    if any(not isinstance(c,(int,F)) for c in h.values()) or not hermitian(h):
        raise ValueError('Obstruction implementation requires real rational Hermitian input')
    if quotient_number_ideal:h=number_ideal_core(h,m)
    g,den,rows,nnz=integer_gram(h,m);gram_seconds=time.monotonic()-start
    # Coordinates Eii,Eij+Eji have domain Frobenius metric diag(1,2).
    metric=[1 if i==j else 2 for i,j in generators(m)]
    shifted=[[x*threshold.denominator for x in row] for row in g]
    for i,d in enumerate(metric):shifted[i][i]-=threshold.numerator*den*den*d
    inertia=inertia_integer(shifted);k=inertia['negative'];extra={}
    if complex_unitary:
        ag,ad,arows,annz=integer_gram(h,m,True)
        ashift=[[x*threshold.denominator for x in row] for row in ag]
        for i in range(len(ag)):ashift[i][i]-=2*threshold.numerator*ad*ad
        ai=inertia_integer(ashift);k+=ai['negative']
        extra={'imaginary_antisymmetric_inertia':ai,'imaginary_unknowns':len(ag),
               'imaginary_equation_rows':arows,'imaginary_map_nonzeros':annz}
    squared=F(max(0,m-k))*threshold/4
    # A downward-rounded rational norm bound avoids relying on floating sqrt.
    precision=10**12;lower=F(isqrt((squared.numerator*precision**2)//squared.denominator),precision)
    return {'modes':m,'unknowns':len(g),'equation_rows':rows,'map_nonzeros':nnz,'threshold':str(threshold),
            'inertia':inertia,'total_negative_inertia':k,'complex_unitary':complex_unitary,
            'quotient_one_body_number_ideal':quotient_number_ideal,
            'frobenius_distance_squared_lower':str(squared),**extra,
            'frobenius_distance_lower':str(lower),'frobenius_distance_lower_float':float(lower),
            'excludes_coefficient_l1_residual_0_0016':lower>F(1,625),
            'gram_seconds':gram_seconds,'total_seconds':time.monotonic()-start,
            'many_body_states_enumerated':0,
            'scope':'quartic coefficient Frobenius distance to any '+('complex unitarily' if complex_unitary else 'real orthogonally')+' rotated density-density interaction'+(' plus arbitrary one-body-multiplier number ideal' if quotient_number_ideal else '')}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threshold',type=F,default=F(1,10000));p.add_argument('--complex-unitary',action='store_true')
    p.add_argument('--quotient-number-ideal',action='store_true');a=p.parse_args()
    f=json.loads(a.fixture.read_text());h=decode(f['hamiltonian'],f['modes'],4)
    r=obstruction(h,f['modes'],a.threshold,a.complex_unitary,a.quotient_number_ideal);r['fixture_sha256']=hashlib.sha256(a.fixture.read_bytes()).hexdigest()
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
