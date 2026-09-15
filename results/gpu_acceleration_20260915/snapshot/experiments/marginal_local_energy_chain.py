"""Global Hubbard energy intervals from recomputed all-Fock local positivity.

The local cluster is K_L(U)=H_L(U)-U/2*(Nhat-L). On the global half-filled
sector the corresponding centered chain is exactly the physical H_N.
"""
from fractions import Fraction as F
import json
from experiments.marginal_local_hubbard_block import replay as replay_block


def _replay(certificate, verify_block):
    if type(certificate) is not dict or certificate.get('kind') != 'local_energy_chain_v1':
        raise ValueError('Unsupported chain certificate')
    N, L = certificate.get('sites'), certificate.get('block_length')
    if type(N) is not int or not 2 <= N <= 10**9 or N % 2:
        raise ValueError('Even chain length in2..10^9 required')
    if type(L) is not int or L not in (2,4,6) or N < L:
        raise ValueError('Local block length must be2,4,6 and no larger than chain')
    for key in ('U','t'):
        if type(certificate.get(key)) not in (int,str,F):
            raise ValueError('Exact chain parameters required')
    U, t = F(certificate['U']), F(certificate['t'])
    if any(x < 0 or x > 10**6 or x.denominator > 10**6 for x in (U,t)):
        raise ValueError('Invalid bounded chain parameters')
    q, r = divmod(N,L)
    required = {L} | ({r} if r else set())
    blocks = certificate.get('blocks')
    if type(blocks) is not list or len(blocks) != len(required):
        raise ValueError('Exactly the required full and remainder block certificates are needed')
    receipts = {}
    for block in blocks:
        if type(block) is not dict or block.get('sites') not in required or block['sites'] in receipts:
            raise ValueError('Missing, duplicate, or unrelated local block')
        if F(block.get('U')) != U or F(block.get('t')) != t:
            raise ValueError('Local model does not match the chain')
        result = verify_block(block)
        if (list(map(F,result['onsite_profile'])) != [U]*block['sites']
                or list(map(F,result['hopping_profile'])) != [t]*(block['sites']-1)):
            raise ValueError('Tiled physical blocks must have uniform chain coefficients')
        if result.get('upper_quotient') is None or F(result['upper_norm']) <= 0:
            raise ValueError('Every tiled block needs a physical positive-norm upper witness')
        receipts[block['sites']] = result
    if set(receipts) != required:
        raise ValueError('Incomplete tiling')
    cuts = q-1+bool(r)
    lower = q*F(receipts[L]['lower'])-2*t*cuts
    upper = q*F(receipts[L]['upper_quotient'])
    if r:
        lower += F(receipts[r]['lower'])
        upper += F(receipts[r]['upper_quotient'])
    tiled_lower = lower
    overlap_lower = None
    overlap_receipt = None
    if 'overlap' in certificate:
        block = certificate['overlap']
        W = block.get('sites') if type(block) is dict else None
        if type(W) is not int or W not in (2,4,6) or N <= W or N < 4:
            raise ValueError('Overlap construction requires window2,4,6 and N>window')
        v = U*F(W-1,W)
        if (type(block) is not dict
                or F(block.get('U')) != v or F(block.get('t')) != t):
            raise ValueError('Overlap window must use onsite U*(L-1)/L and the same hopping')
        overlap_receipt = verify_block(block)
        overlap_lower = F(N,W-1)*F(overlap_receipt['lower'])-2*t
        lower = max(lower,overlap_lower)
    if lower > upper:
        raise ValueError('Global lower exceeds physical product upper')
    return {'accepted':True,'sites':N,'block_length':L,'full_blocks':q,
            'remainder_sites':r,'removed_tiling_bonds':cuts,
            'lower':str(lower),'upper':str(upper),'width':str(upper-lower),
            'lower_per_site':str(lower/N),'upper_per_site':str(upper/N),
            'width_per_site':str((upper-lower)/N),
            'tiling_lower':str(tiled_lower),
            'overlap_lower':str(overlap_lower) if overlap_lower is not None else None,
            'local_replays':{str(k):v for k,v in receipts.items()},
            'overlap_replay':overlap_receipt,
            'scope':'Exact local all-Fock positivity plus analytic tiling/overlap operator identities and fixed-charge product upper states for the uniform open half-filled Hubbard chain. Whole-chain state enumeration is absent; local work still grows exponentially with block length. No general molecular or physical-marginal representability claim.'}


def replay(certificate):
    return _replay(certificate,replay_block)


def replay_many(certificates):
    """Verify repeated local certificates once within this fresh batch only."""
    if type(certificates) is not list or not 1 <= len(certificates) <= 32:
        raise ValueError('One to32 chain certificates required')
    cache = {}
    def verify(block):
        key = json.dumps(block,sort_keys=True,default=str)
        if key not in cache:
            cache[key] = replay_block(block)
        return cache[key]
    results = [_replay(c,verify) for c in certificates]
    return {'accepted':True,'unique_local_certificates':len(cache),'chains':results,
            'cache_scope':'Fresh in-process exact replay cache; no submitted receipts or persisted validation cache accepted.'}
