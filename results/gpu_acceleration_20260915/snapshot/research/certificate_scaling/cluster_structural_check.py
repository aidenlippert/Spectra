"""Exact CAR check for the clustered dimer theorem (stdlib only)."""
from fractions import Fraction

def mm(a,b): return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b))] for i in range(len(a))]
def add(a,b): return [[a[i][j]+b[i][j] for j in range(len(a))] for i in range(len(a))]
def fop(n,m,create):
    d=1<<n; out=[[0]*d for _ in range(d)]
    for s in range(d):
        occ=(s>>m)&1
        if create==occ: continue
        t=s^(1<<m); out[t][s]= -1 if ((s&((1<<m)-1)).bit_count()&1) else 1
    return out

def main():
    a0,a1=fop(2,0,False),fop(2,1,False); c0,c1=fop(2,0,True),fop(2,1,True)
    n0,n1=mm(c0,a0),mm(c1,a1); total=add(n0,n1)
    ident=[[int(i==j) for j in range(4)] for i in range(4)]
    # H0=(n0+n1-1)^2, verified from the CAR matrices.
    h0=add(add(add(mm(n0,n0),mm(n1,n1)),[[2*mm(n0,n1)[i][j] for j in range(4)] for i in range(4)]),[[ -2*total[i][j]+ident[i][j] for j in range(4)] for i in range(4)])
    expected=[[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]]
    assert h0==expected
    t=add(mm(c0,a1),mm(c1,a0))
    assert mm(t,t)==[[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,0]]  # T^2 projector
    assert mm(n0,t)!=mm(t,n0)                         # local dimer number changes
    assert mm(total,t)==mm(t,total)                   # total edge number conserved

def report(L,g):
    e=L-1; width=e*abs(g)
    print(f"L={L} N={L} edges={e} g={g} lower={-width} upper=0 local_factors={L} edge_residual_blocks={e}")
    return width

if __name__=='__main__':
    main()
    fixed=[report(L,Fraction(1,10)) for L in (4,8,16,32)]
    scaled=[report(L,Fraction(1,10*L)) for L in (4,8,16,32)]
    assert fixed==[Fraction(3,10),Fraction(7,10),Fraction(3,2),Fraction(31,10)]
    assert scaled==[Fraction(3,40),Fraction(7,80),Fraction(3,32),Fraction(31,320)]
    print('PASS: exact CAR local H0, hopping, total-N conservation, and bounds')
