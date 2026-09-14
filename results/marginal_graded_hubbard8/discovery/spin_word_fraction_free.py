"""Bounded fraction-free exact solve for the 198-row discovery family."""
from fractions import Fraction as F
from math import lcm


def solve(matrix,rhs):
    m=len(rhs);n=len(matrix[0]) if matrix else 0
    if not 1<=n<=m<=198 or len(matrix)!=m or any(len(row)!=n for row in matrix):
        raise ValueError('Bounded rectangular198-row system required')
    matrix=[list(map(F,row)) for row in matrix];rhs=list(map(F,rhs))
    scales=[lcm(*(matrix[i][j].denominator for i in range(m))) for j in range(n)]
    q=lcm(*(v.denominator for v in rhs))
    a=[[int(matrix[i][j]*scales[j]) for j in range(n)]+[int(rhs[i]*q)] for i in range(m)]
    if any(abs(v).bit_length()>4096 for row in a for v in row):
        raise ValueError('Cleared integer input exceeds4096-bit budget')
    previous=1
    for k in range(n):
        pivot_row=next((i for i in range(k,m) if a[i][k]),None)
        if pivot_row is None:raise ValueError('Dependent proposed basis columns')
        a[k],a[pivot_row]=a[pivot_row],a[k];pivot=a[k][k]
        if abs(pivot).bit_length()>50000:raise ValueError('Intermediate integer budget exceeded')
        for i in range(m):
            if i==k:continue
            factor=a[i][k]
            for j in range(k+1,n+1):
                value,remainder=divmod(pivot*a[i][j]-factor*a[k][j],previous)
                if remainder:raise ArithmeticError('Nonexact fraction-free division')
                a[i][j]=value
            a[i][k]=0
        previous=pivot
    if any(row[-1] for row in a[n:]):raise ValueError('Exact proposed basis is inconsistent')
    weights=[F(a[i][-1]*scales[i],previous*q) for i in range(n)]
    if any(sum(v*w for v,w in zip(row,weights))!=b for row,b in zip(matrix,rhs)):
        raise ValueError('Exact original rational constraints failed')
    return weights
