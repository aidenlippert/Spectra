"""Exact spin-reflection-positive Hubbard energy certificate.

The actual bound uses 70-dimensional rational PSD checks. The dense candidate
has 4,900 amplitudes: this is not an enumeration-free or polynomial-size solver.
Only standard-library arithmetic is used by the accepting checker.
"""
from fractions import Fraction as F
from itertools import combinations
from math import lcm


def one_spin_model(rungs=4):
    if type(rungs) is not int or rungs not in (2,4):
        raise ValueError('Only four-site regression and eight-site target supported')
    sites=2*rungs;particles=rungs
    labels=sorted(sum(1<<i for i in choice) for choice in combinations(range(sites),particles))
    index={x:i for i,x in enumerate(labels)}
    edges=[(2*i,2*i+1) for i in range(rungs)]
    edges += [(2*i+s,2*(i+1)+s) for i in range(rungs-1) for s in (0,1)]
    eta=[(-1)**(i//2+i%2) for i in range(sites)]
    if any(eta[i]*eta[j]!=-1 for i,j in edges):
        raise AssertionError('Particle-hole sign rule requires bipartite edges')
    K=[[0 for _ in labels] for _ in labels]
    for col,label in enumerate(labels):
        for a,b in edges:
            for dst,src in ((a,b),(b,a)):
                if label&(1<<src) and not label&(1<<dst):
                    lo,hi=sorted((dst,src))
                    between=((1<<hi)-1)^((1<<(lo+1))-1)
                    sign=-1 if (label&between).bit_count()%2 else 1
                    moved=label^(1<<dst)^(1<<src)
                    K[index[moved]][col]-=sign
    if any(K[i][j]!=K[j][i] for i in range(len(K)) for j in range(i)):
        raise AssertionError('Single-spin hopping is not Hermitian')
    return labels,K,edges,eta


def apply_map(C,labels,K,U=8):
    """L(C)=K C+C K - U sum_i n_i C n_i+U*Nup*C."""
    n=len(labels);particles=labels[0].bit_count()
    rows=[[(j,x) for j,x in enumerate(row) if x] for row in K]
    return [[sum(x*C[k][j] for k,x in rows[i])
             +sum(x*C[i][k] for k,x in rows[j])
             +U*(particles-(labels[i]&labels[j]).bit_count())*C[i][j]
             for j in range(n)] for i in range(n)]


def integer_psd(matrix,strict=False):
    """Exact fraction-free LDL with null-pivot and symmetry checks."""
    n=len(matrix)
    if not 1<=n<=70 or any(len(row)!=n for row in matrix):
        raise ValueError('PSD dimension outside 1..70')
    if any(type(x) is not int for row in matrix for x in row):
        raise ValueError('Integer PSD entries required')
    if any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(i)):
        raise ValueError('PSD matrix is not symmetric')
    a=[row[:] for row in matrix];previous=1;rank=0;max_bits=0
    for k in range(n):
        pivot=a[k][k]
        if pivot<0 or (strict and pivot==0):
            raise ValueError('Nonpositive strict pivot or negative PSD pivot')
        if pivot==0:
            if any(a[k][j] for j in range(k+1,n)):
                raise ValueError('Nonzero coupling from a null PSD pivot')
            continue
        max_bits=max(max_bits,pivot.bit_length())
        for i in range(k+1,n):
            for j in range(i,n):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                if numerator%previous:
                    raise ValueError('Nonexact fraction-free division')
                a[i][j]=a[j][i]=numerator//previous
        previous=pivot;rank+=1
    return {'dimension':n,'rank':rank,'strict':strict,'max_pivot_bits':max_bits}


def verify(payload):
    if (not isinstance(payload,dict) or payload.get('kind')!='half_filled_bipartite_hubbard_positive_cone_v1'
            or payload.get('rungs')!=4 or payload.get('U')!='8' or payload.get('t')!='1'
            or payload.get('target_spin_populations')!=[4,4]):
        raise ValueError('Certificate does not bind the fixed physical model')
    labels,K,edges,eta=one_spin_model()
    if payload.get('one_spin_labels')!=labels:
        raise ValueError('Wrong one-spin ordering')
    C=payload['integer_C'];n=len(labels)
    if not isinstance(C,list) or len(C)!=n or any(not isinstance(row,list) or len(row)!=n for row in C):
        raise ValueError('Wrong amplitude matrix shape')
    c_check=integer_psd(C,strict=True)
    if type(payload['lower_over_t']) not in (str,int):
        raise ValueError('Exact rational lower required')
    lower=F(payload['lower_over_t']);p,q=lower.numerator,lower.denominator
    action=apply_map(C,labels,K)
    residual=[[q*action[i][j]-p*C[i][j] for j in range(n)] for i in range(n)]
    r_check=integer_psd(residual)
    norm=sum(x*x for row in C for x in row)
    numerator=sum(C[i][j]*action[i][j] for i in range(n) for j in range(n))
    upper=F(numerator,norm)
    if lower>upper:raise AssertionError('Certified endpoints inconsistent')
    width=upper-lower
    return {'status':'accepted_exact_positive_cone_interval',
            'energy_lower_over_t':str(lower),'energy_upper_over_t':str(upper),
            'width_over_t':str(width),'width_float':float(width),
            'target_width_over_t':'1/1000','meets_target':width<=F(1,1000),
            'amplitude_positive_definite':c_check,'cone_residual_psd':r_check,
            'one_spin_configurations_enumerated':n,'amplitude_entries':n*n,
            'full_spinful_Hamiltonian_constructed':False,
            'dense_spinful_amplitude_information_retained':True,
            'particle_hole_constant_over_t':32,'physical_hopping_edges':edges,
            'scope':'Fixed eight-site half-filled repulsive bipartite Hubbard; also full fixed-N=8 spin ground by SU(2); not a compact/scalable general solver'}


def original_to_matrix_label(state,rungs=4):
    """Exact signed CAR particle-hole basis mapping, for regression only."""
    sites=2*rungs
    up=[i for i in range(sites) if state&(1<<(2*i))]
    down=[i for i in range(sites) if state&(1<<(2*i+1))]
    # Interleaved creation order -> all up creations, then all down creations.
    crossings=sum(j>i for i in down for j in up)
    eta=[(-1)**(i//2+i%2) for i in range(sites)]
    phase=(-1)**(crossings+sum(down))
    for i in down:phase*=eta[i]
    a=sum(1<<i for i in up)
    b=((1<<sites)-1)^sum(1<<i for i in down)
    return a,b,phase
