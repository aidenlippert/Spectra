"""Exact GF(2) parity masks commuting with a real CAR Hamiltonian."""
from fractions import Fraction


from experiments.marginal_symbolic import canonical,validate_word


def spin_parity(poly, modes):
    if type(modes) is not int or modes <= 0 or modes % 2: raise ValueError("modes must be positive even")
    if not isinstance(poly, dict): raise ValueError("polynomial must be dict")
    rows=[]; half=modes//2
    for w,c in poly.items():
        if type(c) is bool or not isinstance(c,(int,Fraction)): raise ValueError("coefficients must be exact real")
        validate_word(w,modes,4)
        if type(w) is not tuple or canonical({w:Fraction(1)})!={w:Fraction(1)}:raise ValueError('Noncanonical word')
        if not c:continue
        row=0
        for _,i in w: row ^= 1 << (i//2)
        if row: rows.append(row)
    # nullspace over GF(2), with pivot elimination and free-variable basis.
    piv=[]
    for col in range(half):
        p=next((r for r in range(len(piv),len(rows)) if (rows[r]>>col)&1),None)
        if p is None: continue
        rows[len(piv)],rows[p]=rows[p],rows[len(piv)]
        for r in range(len(rows)):
            if r!=len(piv) and ((rows[r]>>col)&1): rows[r]^=rows[len(piv)]
        piv.append(col)
    free=[c for c in range(half) if c not in piv]; basis=[]
    for f in free:
        x=1<<f
        for r,col in reversed(list(enumerate(piv))):
            if (rows[r] & x).bit_count() & 1: x |= 1<<col
        mask=0
        for s in range(half):
            if x>>s&1: mask |= (1<<(2*s))|(1<<(2*s+1))
        basis.append(mask)
    return tuple(basis)
