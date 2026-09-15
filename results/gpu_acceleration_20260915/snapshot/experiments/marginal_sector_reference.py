"""Exact ground-energy brackets for the matched, half-filled hopping model.

All conserved pair-occupation sectors are covered. Rational leading-minor
recurrences certify the smallest eigenvalue of each collective Jacobi matrix.
No numerical libraries or occupation-basis enumeration are used here.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_symbolic import decode


def jacobi_data(modes, t, doubles):
    if type(modes) is not int or modes < 4 or modes % 2:
        raise ValueError("Even M>=4 required")
    half = modes//2
    if type(doubles) is not int or not 0 <= doubles <= half//2:
        raise ValueError("Invalid number of doubly occupied pairs")
    t = F(t)
    singles = half-2*doubles
    diagonal = [F((doubles+k)*(doubles+k-1)//2
                  +(doubles+singles-k)*(doubles+singles-k-1)//2)
                for k in range(singles+1)]
    offdiagonal_squared = [t*t*(k+1)*(singles-k) for k in range(singles)]
    return diagonal, offdiagonal_squared


def above_lower_bound(diagonal, offdiagonal_squared, bound):
    """True iff the real symmetric Jacobi matrix minus bound*I is PD."""
    if not diagonal or len(offdiagonal_squared) != len(diagonal)-1:
        raise ValueError("Inconsistent Jacobi data")
    if any(q < 0 for q in offdiagonal_squared):
        raise ValueError("Squared offdiagonals must be nonnegative")
    previous = F(1)
    current = F(diagonal[0])-bound
    if current <= 0:
        return False
    for k in range(1,len(diagonal)):
        following = (diagonal[k]-bound)*current-offdiagonal_squared[k-1]*previous
        if following <= 0:
            return False
        previous, current = current, following
    return True


def bracket(modes, t, doubles, tolerance=F(1,10**12)):
    tolerance = F(tolerance); t = F(t)
    if tolerance <= 0:
        raise ValueError("Positive tolerance required")
    diagonal, squared = jacobi_data(modes,t,doubles)
    singles = len(diagonal)-1
    # The hopping hypercube has operator norm |t|*singles. A diagonal
    # basis vector supplies the upper bound min(diagonal).
    lower = min(diagonal)-abs(t)*singles-1
    upper = min(diagonal)
    assert above_lower_bound(diagonal,squared,lower)
    assert not above_lower_bound(diagonal,squared,upper)
    while upper-lower > tolerance:
        middle = (lower+upper)/2
        if above_lower_bound(diagonal,squared,middle):
            lower = middle
        else:
            upper = middle
    return {"doubles": doubles, "singles": singles, "dimension": singles+1,
            "lower": str(lower), "upper": str(upper),
            "lower_float": float(lower), "upper_float": float(upper)}


def reference(modes, t, tolerance=F(1,10**12)):
    # Validate even mode count before deriving sector range.
    jacobi_data(modes,t,0)
    sectors = [bracket(modes,t,d,tolerance) for d in range(modes//4+1)]
    lower = min(F(s['lower']) for s in sectors)
    upper = min(F(s['upper']) for s in sectors)
    winner = min(sectors,key=lambda s:F(s['upper']))
    isolated = all(F(winner['upper']) < F(s['lower']) for s in sectors if s is not winner)
    return {"modes": modes, "particles": modes//2, "t": str(F(t)), "sectors": sectors,
            "lower": str(lower), "upper": str(upper), "width": str(upper-lower),
            "lower_float": float(lower), "upper_float": float(upper),
            "isolated_ground_sector_doubles": winner['doubles'] if isolated else None}


def for_certificate(certificate, tolerance=F(1,10**12)):
    modes = certificate['modes']; particles = certificate['particles']
    t = F(certificate['variational_upper']['t'])
    if particles != modes//2 or decode(certificate['hamiltonian'],modes,4) != hopping_polynomial(modes,t):
        raise ValueError("Reference requires the exact matched half-filled Hamiltonian")
    return reference(modes,t,tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes',type=int,default=10)
    parser.add_argument('--t',default='1/5')
    parser.add_argument('--certificate',type=Path)
    args = parser.parse_args()
    result = (for_certificate(json.loads(args.certificate.read_text())) if args.certificate
              else reference(args.modes,F(args.t)))
    print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
