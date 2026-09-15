"""Rational l1 disks for bounded complex exponential-polynomial readout.

Exact coefficient aggregation precedes shared exponential evaluation. Scaling,
Taylor truncation and every grid rounding are enclosed; unsupported precision
or resource requests fail rather than returning an unproved bound.
"""
from fractions import Fraction as F
from .v7_certificate import rational, BudgetExceeded

MAX_TERMS=50000
MAX_DEGREE=32
MAX_TAYLOR=128
MAX_SCALINGS=16


def plus(a,b): return rational(a[0]+b[0]),rational(a[1]+b[1])
def times(a,b): return rational(a[0]*b[0]-a[1]*b[1]),rational(a[0]*b[1]+a[1]*b[0])
def scale(a,b): return rational(a[0]*b),rational(a[1]*b)
def norm(a): return rational(abs(a[0])+abs(a[1]))


def rounded_disk(mid,radius,bits,cost):
    grid=1<<bits
    out=tuple(F((v*grid).numerator//(v*grid).denominator,grid) for v in mid)
    radius=rational(radius+norm((mid[0]-out[0],mid[1]-out[1])))
    z=radius*grid
    radius=F(-(-z.numerator//z.denominator),grid)
    cost['grid_roundings']+=3
    return out,rational(radius)


def exp_disk(z,tolerance,cost):
    """Return exp(z) with l1 error <= tolerance, or refuse."""
    mag=norm(z)
    if not mag:return (F(1),F(0)),F(0)
    scalings=0
    while mag/F(1<<scalings)>F(1,2):
        scalings+=1
        if scalings>MAX_SCALINGS:raise BudgetExceeded('exponential scaling budget')
    # Verify the final disk rather than assuming small errors stay small after
    # squaring. Every grid change is charged to the outward l1 error radius.
    bits=max(64,tolerance.denominator.bit_length()-abs(tolerance.numerator).bit_length()+24+8*scalings)
    if bits>4096:raise BudgetExceeded('readout grid precision budget')
    w=scale(z,F(1,1<<scalings));wm=norm(w)
    target=rational(tolerance/F(16*(1<<scalings)))
    for refinement in range(16):
        term=(F(1),F(0));mid=term
        for n in range(1,MAX_TAYLOR+1):
            term=scale(times(term,w),F(1,n));mid=plus(mid,term)
            ratio=wm/F(n+1)
            radius=rational(norm(term)*ratio/(1-ratio))
            cost['taylor_terms']+=1
            if radius<=target:break
        else:raise BudgetExceeded('readout Taylor iteration budget')
        mid,radius=rounded_disk(mid,radius,bits,cost)
        for _ in range(scalings):
            radius=rational(2*norm(mid)*radius+radius**2)
            mid=times(mid,mid)
            mid,radius=rounded_disk(mid,radius,bits,cost)
            cost['squarings']+=1
        if radius<=tolerance:return mid,radius
        cost['precision_refinements']+=1
        target=rational(target/16)
        bits+=16
        if bits>4096:raise BudgetExceeded('readout grid precision budget')
    raise BudgetExceeded('readout refinement budget')


def evaluate_mode_map(mapping,time,tolerance):
    """Return word->midpoint pairs, SUM of word l1 errors, and work counts.

    Rows: (real exponent, imaginary exponent, degree, real/imag coefficient).
    No float/string coercion or unbounded list creation.
    """
    time,tolerance=rational(time),rational(tolerance)
    if time<0 or tolerance<=0:raise ValueError('readout time/tolerance')
    if not isinstance(mapping,dict) or len(mapping)>512:raise BudgetExceeded('readout support budget')
    cost=dict(input_terms=0,unique_exponents=0,taylor_terms=0,squarings=0,
              precision_refinements=0,grid_roundings=0,coefficient_aggregations=0)
    aggregated={};powers={0:F(1)}
    for word,rows in mapping.items():
        if not isinstance(word,str):raise ValueError('readout key must be a string')
        for row in rows:
            cost['input_terms']+=1
            if cost['input_terms']>MAX_TERMS:raise BudgetExceeded('readout entry budget')
            if not isinstance(row,(tuple,list)) or len(row)!=5:raise ValueError('readout mode fields')
            a,b,k,cr,ci=row
            a,b,cr,ci=map(rational,(a,b,cr,ci))
            if a>0 or type(k)is not int or not 0<=k<=MAX_DEGREE:raise ValueError('readout exponent/degree')
            if k not in powers:powers[k]=rational(time**k)
            key=(a,b);group=aggregated.setdefault(key,{})
            c=plus(group.get(word,(F(0),F(0))),scale((cr,ci),powers[k]))
            if c!=(F(0),F(0)):group[word]=c
            elif word in group:del group[word]
            cost['coefficient_aggregations']+=1
    aggregated={key:g for key,g in aggregated.items() if g}
    cost['unique_exponents']=len(aggregated)
    result={word:(F(0),F(0)) for word in mapping};error=F(0)
    for (a,b),group in aggregated.items():
        weight=rational(sum((norm(c) for c in group.values()),F(0)))
        target=rational(tolerance/F(max(1,len(aggregated)))/max(F(1),weight))
        mid,radius=exp_disk((rational(a*time),rational(b*time)),target,cost)
        for word,c in group.items():result[word]=plus(result[word],times(c,mid))
        error=rational(error+weight*radius)
    if error>tolerance:raise ValueError('internal readout error allocation failure')
    return result,error,cost


def evaluate_modes(terms,time,tolerance):
    result,error,cost=evaluate_mode_map({'value':terms},time,tolerance)
    return result['value'],error,cost
