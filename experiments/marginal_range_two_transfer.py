"""Exact physical filtered-state transfer with next-nearest density terms.

Only the energy tensors change. Each range-two term touches at most one cut
filter for blocks of at least four sites, and at most two adjacent blocks.
The existing norm and 48-dimensional energy recurrence therefore applies.
"""
from fractions import Fraction as F

from experiments.marginal_local_hubbard_block import _actions, _exact
from experiments.marginal_range_two_density import diagonal_value
from experiments.marginal_boundary_transfer import (
    compile_state as compile_nearest, contract as contract_nearest, enclose,
    _partial, _edge_tensor, _compose, _add,
)
from experiments.marginal_boundary_unitary import _physical_source


def _dressed_range_two(coupling):
    h = _actions(2, 0, 1)
    identity = {s: {s: 1} for s in range(16)}
    powers = [identity, {s: {t: -a for t, a in image.items()} for s, image in h.items()}, _compose(h, h)]
    # Six-site patch: the filter acts on sites2 and3. Four range-two pairs
    # touch one of these sites: (0,2),(1,3),(2,4),(3,5).
    embedded = [{s: {(s & ~240)+(t << 4): a for t, a in op[(s >> 4) & 15].items()}
                 for s in range(4096)} for op in powers]
    diagonal = {s: {s: diagonal_value(s, [coupling]*4)} for s in range(4096)}
    for m in range(3):
        for n in range(m, 3):
            energy = _compose(embedded[m], _compose(diagonal, embedded[n]))
            if m != n:
                energy = _add(energy, _compose(embedded[n], _compose(diagonal, embedded[m])))
            if any(energy[t].get(s, 0) != a for s, image in energy.items() for t, a in image.items()):
                raise ValueError('Dressed range-two patch must be exactly Hermitian')
            yield m, n, energy


def compile_state(state, sites, U=4, t=1, V=0, W=0):
    """Compile from an explicit physical integer state, never supplied RDMs."""
    W = _exact(W)
    compiled = compile_nearest(state, sites, U, t, V)
    state = {s: a for s, a in state.items() if a}
    compiled['target']['W'] = str(W)
    if not W:
        return compiled
    weights = {
        'J': [0]+[W]*(sites-4)+[0],
        'L': [W]+[0]*(sites-3),
        'R': [0]*(sites-3)+[W],
    }
    for name, profile in weights.items():
        image = {s: a*diagonal_value(s, profile) for s, a in state.items()}
        addition = _edge_tensor(image, state, sites)
        compiled[name] = [[a+addition[i][j] for j, a in enumerate(row)]
                          for i, row in enumerate(compiled[name])]
    # Paired exterior boundary indices remain the16 transfer coordinates.
    rho_a = _partial(state, state, (0, sites-3, sites-2, sites-1), sites)
    rho_b = _partial(state, state, (0, 1, 2, sites-1), sites)
    left, right = {}, {}
    for (s, bra), value in rho_a.items():
        left.setdefault((s >> 2, bra >> 2), []).append((s % 4+4*(bra % 4), value))
    for (s, bra), value in rho_b.items():
        right.setdefault((s % 64, bra % 64), []).append((s // 64+4*(bra // 64), value))
    kernels = []
    for old, (m, n, energy) in zip(compiled['kernels'], _dressed_range_two(W)):
        if old[:2] != (m, n):
            raise ValueError('Cut-polynomial kernel order differs')
        addition = [[F(0)]*16 for _ in range(16)]
        for s, image in energy.items():
            for bra, coefficient in image.items():
                for i, x in left.get((s % 64, bra % 64), []):
                    for j, y in right.get((s // 64, bra // 64), []):
                        addition[i][j] += coefficient*x*y
        kernels.append((m, n, old[2], [[a+addition[i][j] for j, a in enumerate(row)]
                                     for i, row in enumerate(old[3])]))
    if len(kernels) != 6:
        raise ValueError('Incomplete range-two cut polynomial')
    compiled['kernels'] = kernels
    compiled['range_two_patch_sites'] = 6
    return compiled


def compile_block(c, h, U=4, t=1, V=0, W=0):
    oracle, state, checked = _physical_source(c, h)
    full = {s: a*phase for r, a in state.items() for s, phase in oracle.orbit(r)[3].items()}
    result = compile_state(full, 8, U, t, V, W)
    result['source_upper_replay'] = checked
    return result


def contract(compiled, a, b, blocks):
    result = contract_nearest(compiled, a, b, blocks)
    result['scope'] = 'Exact finite physical norm and U,t,V,W energy, including next-nearest density. Each dressed term occupies at most two adjacent physical blocks. The existing recurrence remains capped at64 blocks; no generic molecular transfer claim.'
    return result
