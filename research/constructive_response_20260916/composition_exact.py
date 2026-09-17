"""Small exact matrix regressions for established response/composition rules.

These finite identities do not establish compactness or a complexity theorem.
"""
from fractions import Fraction as F


def shape(a):
    if not a or not a[0] or any(len(r) != len(a[0]) for r in a):
        raise ValueError('Nonempty rectangular matrix required')
    if any(not isinstance(x,(int,F)) for row in a for x in row):
        raise ValueError('Exact rational entries required')
    return len(a),len(a[0])


def mm(a,b):
    na,ma=shape(a);nb,mb=shape(b)
    if ma != nb:
        raise ValueError('Incompatible matrix product')
    return [[sum((a[i][k]*b[k][j] for k in range(ma)),F(0))
             for j in range(mb)] for i in range(na)]


def tr(a):
    n,m=shape(a)
    return [[a[j][i] for j in range(n)] for i in range(m)]


def add(a,b):
    if shape(a) != shape(b):
        raise ValueError('Incompatible sum')
    return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def sc(c,a):
    shape(a)
    return [[F(c)*x for x in row] for row in a]


def sub(a,b):
    return add(a,sc(-1,b))


def eye(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]


def is_psd(a,strict=False):
    n,m=shape(a)
    if n != m or a != tr(a):
        return False
    residual=[[F(x) for x in row] for row in a]
    for k in range(n):
        pivot=residual[k][k]
        if pivot < 0 or (strict and pivot == 0):
            return False
        if pivot == 0:
            if any(residual[j][k] for j in range(k+1,n)):
                return False
            continue
        for i in range(k+1,n):
            for j in range(i,n):
                value=residual[i][j]-residual[i][k]*residual[j][k]/pivot
                residual[i][j]=residual[j][i]=value
    return True


def inv(a):
    n,m=shape(a)
    if n != m:
        raise ValueError('Square inverse required')
    z=[[F(x) for x in row]+eye(n)[i] for i,row in enumerate(a)]
    for j in range(n):
        candidates=[k for k in range(j,n) if z[k][j]]
        if not candidates:
            raise ValueError('Singular matrix')
        k=candidates[0];z[j],z[k]=z[k],z[j]
        q=z[j][j];z[j]=[x/q for x in z[j]]
        for i in range(n):
            if i != j:
                q=z[i][j];z[i]=[x-q*y for x,y in zip(z[i],z[j])]
    result=[row[n:] for row in z]
    if mm(a,result) != eye(n):
        raise AssertionError('Inverse identity failed')
    return result


def block(a,rows,cols):
    return [[a[i][j] for j in cols] for i in rows]


def schur(a,retained):
    n,m=shape(a)
    if n != m or a != tr(a):
        raise ValueError('Symmetric square operator required')
    retained=list(retained)
    if not retained or len(set(retained)) != len(retained) or any(i<0 or i>=n for i in retained):
        raise ValueError('Invalid retained indices')
    eliminated=[i for i in range(n) if i not in retained]
    if not eliminated:
        return block(a,retained,retained)
    C=block(a,eliminated,eliminated)
    if not is_psd(C,strict=True):
        raise ValueError('Eliminated block is not positive definite')
    B=block(a,eliminated,retained)
    return sub(block(a,retained,retained),mm(mm(tr(B),inv(C)),B))


def response_bounds(A,B,X,D,eta):
    eta=F(eta)
    if not 0 <= eta < 1 or not is_psd(D,strict=True):
        raise ValueError('Positive preconditioner and eta in [0,1) required')
    if not is_psd(sub(A,sc(1-eta,D))) or not is_psd(sub(sc(1+eta,D),A)):
        raise ValueError('Relative preconditioner sandwich failed')
    R=sub(B,mm(A,X))
    base=sub(add(mm(tr(X),B),mm(tr(B),X)),mm(mm(tr(X),A),X))
    correction=mm(mm(tr(R),inv(D)),R)
    return (add(base,sc(1/(1+eta),correction)),
            add(base,sc(1/(1-eta),correction)))


def exact_composition():
    A=[[F(3),F(1,2)],[F(1,2),F(2)]]
    D=[[F(3),F(0)],[F(0),F(2)]]
    B=[[F(1),F(1,2)],[F(0),F(1)]]
    X=[[F(1,4),F(0)],[F(0),F(1,2)]]
    assert mm(A,D) != mm(D,A)
    L,U=response_bounds(A,B,X,D,F(1,4))
    exact=mm(mm(tr(B),inv(A)),B)
    assert is_psd(sub(exact,L)) and is_psd(sub(U,exact))
    return {'response_exact':exact,'response_lower':L,'response_upper':U}


def nested_exact():
    R=[[F(1),F(1,2),F(-1,3),F(1,4),F(1,5)],
       [F(1,3),F(2),F(1,2),F(-1,4),F(2,5)],
       [F(1,2),F(-1,3),F(1),F(1,2),F(-1,5)],
       [F(1,5),F(1,6),F(1,3),F(2),F(1,2)],
       [F(-1,4),F(1,7),F(1,6),F(1,3),F(1)]]
    M=add(eye(5),mm(tr(R),R))
    v=[[F(1),F(-1),F(2),F(1),F(-2)]]
    error=sc(F(1,50),mm(tr(v),v))
    low=sub(M,error);high=add(M,error)
    assert is_psd(low,strict=True)
    direct=schur(M,[0,1])
    first=schur(M,[0,1,2]);nested=schur(first,[0,1])
    lower=schur(schur(low,[0,1,2]),[0,1])
    upper=schur(schur(high,[0,1,2]),[0,1])
    assert direct == nested
    assert lower == schur(low,[0,1]) and upper == schur(high,[0,1])
    assert is_psd(sub(direct,lower)) and is_psd(sub(upper,direct))
    return {'direct':direct,'nested':nested,'lower':lower,'upper':upper}


def additive_error_counterexample():
    # A unit change in one eliminated diagonal causes an eighteen-unit change
    # in the Schur complement. Both full matrices here are positive definite.
    low=[[F(40),F(6)],[F(6),F(1)]]
    high=add(low,[[F(0),F(0)],[F(0),F(1)]])
    assert is_psd(low,strict=True)
    change=sub(schur(high,[0]),schur(low,[0]))[0][0]
    assert change == 18
    return change
