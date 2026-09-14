"""Bounded physical H8 polynomial witness and extensive product-state uppers."""
from fractions import Fraction as F

from experiments.marginal_local_hubbard_block import _terms, _exact
from experiments.marginal_symbolic import canonical, decode, encode
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
from experiments.marginal_transfer_verify import apply_word
from results.marginal_graded_hubbard8.discovery.singlet_moment_energy import (
    upper_monomials, polynomial_upper)


def _uniform_h(sites):
    # Generate O(sites) CAR words only: no local Fock matrix is constructed.
    return {'modes':2*sites,'particles':sites,
            'hamiltonian':encode(canonical(dict(_terms(sites,F(4),F(1)))))}


def _vector(vector,sites,support):
    if type(vector) is not dict or not 1 <= len(vector) <= support:
        raise ValueError('Bounded nonempty boundary support required')
    parsed = {}
    for key,value in vector.items():
        if type(key) not in (int,str):
            raise ValueError('Integer state label required')
        state = int(key)
        if (state in parsed or not 0 <= state < 4**sites or state.bit_count() != sites
                or type(value) is not int or abs(value) > 10**12):
            raise ValueError('Bounded integer fixed-particle boundary required')
        parsed[state] = value
    if not any(parsed.values()):
        raise ValueError('Nonzero physical boundary required')
    return parsed


def replay_upper8(c,h):
    """Recompute a physical polynomial quotient; no ground-spin theorem needed."""
    if type(c) is not dict or type(c.get('embedding')) is not dict:
        raise ValueError('Certificate embedding required')
    if (type(h) is not dict or type(h.get('modes')) is not int or h['modes'] != 16
            or type(h.get('particles')) is not int or h['particles'] != 8):
        raise ValueError('Original16-mode8-particle Hamiltonian required')
    original = _uniform_h(8)
    if decode(h.get('hamiltonian'),16,4) != decode(original['hamiltonian'],16,4):
        raise ValueError('Hamiltonian must be the uniform open H8 U4 t1 operator')
    vectors = c['embedding'].get('basis')
    table = c.get('upper_chebyshev_coefficients')
    if type(vectors) is not list or not 1 <= len(vectors) <= 14:
        raise ValueError('One to14 physical boundaries required')
    if (type(table) is not list or not 1 <= len(table) <= 12
            or any(type(row) is not list or len(row) != len(vectors) for row in table)):
        raise ValueError('Bounded upper coefficient table shape required')
    if any(type(x) is not int or abs(x) > 10**15 for row in table for x in row):
        raise ValueError('Bounded integer polynomial coefficients required')
    if not any(x for row in table for x in row):
        raise ValueError('Nonzero polynomial recipe required')
    parsed = [_vector(v,8,70) for v in vectors]
    oracle = SymmetryMomentOracle(original)
    # The oracle checks every nonzero input orbit and exact Hamiltonian
    # symmetries, preserving its existing order, source, and support limits.
    moments,work = oracle.moments(parsed,2*len(table)-1)
    upper,norm = polynomial_upper(moments,upper_monomials(table))
    return {'accepted':True,'sites':8,'particles':8,'U':'4','t':'1',
            'upper':str(upper),'norm':str(norm),'moment_work':work,
            'scope':'Fresh physical H8 polynomial Rayleigh quotient only. No supplied moment table, complement-gap proof, or ground-spin theorem is used.'}


def remainder_upper(c,sites):
    if type(sites) is not int or sites not in (2,4,6):
        raise ValueError('Remainder must have2,4,6 sites')
    if type(c) is not dict or type(c.get('sites')) is not int or c['sites'] != sites:
        raise ValueError('Matching remainder certificate required')
    if (_exact(c.get('U')) != 4 or _exact(c.get('t')) != 1
            or c.get('onsite_profile') is not None or c.get('hopping_profile') is not None):
        raise ValueError('Uniform U4 t1 remainder required')
    vector = _vector(c.get('upper_vector'),sites,4**sites)
    norm = sum(a*a for a in vector.values())
    numerator = F(0)
    for source,a in vector.items():
        for word,b in _terms(sites,F(4),F(1)):
            result = apply_word(word,source)
            if result:
                target,phase = result
                numerator += a*b*phase*vector.get(target,0)
    return {'sites':sites,'particles':sites,'norm':str(norm),'upper':str(numerator/norm)}


def replay_tiling(c,h,sites,remainder_cert=None):
    if type(sites) is not int or not 8 <= sites <= 10**9 or sites%2:
        raise ValueError('Even chain size between8 and10^9 required')
    q,r = divmod(sites,8)
    if not r and remainder_cert is not None:
        raise ValueError('No remainder certificate is used for a multiple of8')
    remainder = remainder_upper(remainder_cert,r) if r else None
    base = replay_upper8(c,h)
    total = q*F(base['upper'])+(F(remainder['upper']) if remainder else 0)
    return {'accepted':True,'sites':sites,'blocks':q,'remainder':r,
            'upper':str(total),'upper_per_site':str(total/sites),
            'block_upper':base,'remainder_upper':remainder,
            'interblock_hopping_expectation':'0',
            'scope':'Exact open-chain variational upper from a product of fixed-particle physical blocks. Every hopping term crossing a block boundary changes both block particle numbers and has zero expectation. Global state is specified by the block recipe; it is not expanded.'}
