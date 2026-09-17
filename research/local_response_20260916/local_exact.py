"""Exact local primal/dual certificates, no global determinant construction.

Fixed model: open 2x4 Hubbard ladder, U=8, t=1, Nup=Ndown=4.
Local Fock spaces are enumerated, with all charge sectors included.
"""
from fractions import Fraction as F
from math import lcm


def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def spin_counts(state, sites):
    return (sum((state >> (2*i)) & 1 for i in range(sites)),
            sum((state >> (2*i+1)) & 1 for i in range(sites)))


def sector_labels(sites):
    if sites not in (2, 4):
        raise ValueError('Only local rungs and two-rung patches admitted')
    result = {}
    for state in range(1 << (2*sites)):
        result.setdefault(spin_counts(state, sites), []).append(state)
    return dict(sorted(result.items()))


def hop(state, dst, src):
    if not state & (1 << src) or state & (1 << dst):
        return None
    sign = -1 if (state & ((1 << src)-1)).bit_count() % 2 else 1
    mid = state ^ (1 << src)
    if (mid & ((1 << dst)-1)).bit_count() % 2:
        sign = -sign
    return mid | (1 << dst), sign


def rung_matrix():
    out = zero(16)
    for col in range(16):
        doublons = sum(((col >> (2*i)) & 3) == 3 for i in range(2))
        out[col][col] = F(8*doublons-4*(col.bit_count()-2))
        for spin in (0, 1):
            for dst, src in ((spin, 2+spin), (2+spin, spin)):
                moved = hop(col, dst, src)
                if moved is not None:
                    row, sign = moved
                    out[row][col] -= sign
    return out


def parse_matrix(matrix, n):
    if not isinstance(matrix, list) or len(matrix) != n:
        raise ValueError('Wrong matrix dimension')
    if any(not isinstance(row, list) or len(row) != n for row in matrix):
        raise ValueError('Wrong matrix shape')
    if any(type(x) not in (str, int, F) for row in matrix for x in row):
        raise ValueError('Exact rational entries required')
    result = [list(map(F, row)) for row in matrix]
    if any(result[i][j] != result[j][i] for i in range(n) for j in range(i)):
        raise ValueError('Matrix must be exactly symmetric')
    return result


def check_messages(messages):
    if not isinstance(messages, list) or len(messages) != 2:
        raise ValueError('Exactly two boundary messages required')
    result = [parse_matrix(m, 16) for m in messages]
    counts = [spin_counts(i, 2) for i in range(16)]
    if any(m[i][j] and counts[i] != counts[j]
           for m in result for i in range(16) for j in range(16)):
        raise ValueError('Messages must preserve both spin populations')
    return result


def local_blocks(messages, rungs=4):
    """Three dictionaries: charge key -> (local labels, exact matrix)."""
    if rungs != 4:
        raise ValueError('This model is fixed to four rungs')
    messages = check_messages(messages)
    result = []
    rung = rung_matrix()
    for block_id in range(3):
        blocks = {}
        wl = F(1) if block_id == 0 else F(1, 2)
        wr = F(1) if block_id == 2 else F(1, 2)
        for key, labels in sector_labels(4).items():
            index = {x:i for i, x in enumerate(labels)}
            matrix = zero(len(labels))
            for col, state in enumerate(labels):
                left, right = state & 15, state >> 4
                for new in range(16):
                    vl = wl*rung[new][left]
                    vr = wr*rung[new][right]
                    if block_id:
                        vl -= messages[block_id-1][new][left]
                    if block_id < 2:
                        vr += messages[block_id][new][right]
                    if vl:
                        matrix[index[new + 16*right]][col] += vl
                    if vr:
                        matrix[index[left + 16*new]][col] += vr
                # Exactly two legs and both spins; each direction once.
                for a, b in ((0, 2), (1, 3)):
                    for spin in (0, 1):
                        for dst, src in ((2*a+spin, 2*b+spin),
                                         (2*b+spin, 2*a+spin)):
                            moved = hop(state, dst, src)
                            if moved is not None:
                                target, sign = moved
                                matrix[index[target]][col] -= sign
            if any(matrix[i][j] != matrix[j][i]
                   for i in range(len(labels)) for j in range(i)):
                raise AssertionError('Local CAR Hamiltonian is not symmetric')
            blocks[key] = (labels, matrix)
        result.append(blocks)
    return result


def require_psd(matrix):
    # Fraction-free LDL, the same elimination as marginal_polynomial_sos,
    # kept here without that module's unrelated oracle import chain.
    n = len(matrix)
    if not 1 <= n <= 36 or any(len(row) != n for row in matrix):
        raise ValueError('Local PSD matrix dimension must be in 1..36')
    if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(i)):
        raise ValueError('PSD matrix must be symmetric')
    denominator = lcm(*(x.denominator for row in matrix for x in row))
    a = [[int(x*denominator) for x in row] for row in matrix]
    previous = 1
    for k in range(n):
        pivot = a[k][k]
        if pivot < 0:
            raise ValueError('Negative local PSD pivot')
        if pivot == 0:
            if any(a[k][j] for j in range(k+1, n)):
                raise ValueError('Nonzero coupling at null local PSD pivot')
            continue
        for i in range(k+1, n):
            for j in range(i, n):
                numerator = pivot*a[i][j]-a[i][k]*a[k][j]
                if numerator % previous:
                    raise ValueError('Nonexact fraction-free division')
                a[i][j] = a[j][i] = numerator//previous
        previous = pivot


def require_model(payload):
    if (payload.get('kind') != 'local_hubbard_boundary_v1'
            or payload.get('rungs') != 4 or payload.get('U') != '8'
            or payload.get('t') != '1'):
        raise ValueError('Wrong fixed Hubbard model')


def verify(payload):
    require_model(payload)
    messages = check_messages(payload['messages'])
    shifts = payload['lower_shifts']
    if (not isinstance(shifts, list) or len(shifts) != 3
            or any(type(x) not in (str, int, F) for x in shifts)):
        raise ValueError('Three exact local energy shifts required')
    shifts = list(map(F, shifts))
    checked = 0
    for i, blocks in enumerate(local_blocks(messages)):
        for labels, matrix in blocks.values():
            a = [row[:] for row in matrix]
            for j in range(len(a)):
                a[j][j] -= shifts[i]
            require_psd(a)
            checked += 1
    return {'status':'accepted_exact_local_lower_component',
            'energy_lower_over_t':str(sum(shifts)),
            'local_sector_blocks_checked':checked,
            'local_fock_dimension':256, 'max_psd_dimension':36,
            'global_determinants_enumerated':0,
            'fixed_global_charge':[4,4],
            'chemical_shift':'-4(N-8), zero in the target sector',
            'boundary_identity':'Each boundary message cancels once with its negative'}


def partial_trace(blocks, keep):
    """Blocks map a charge key to (labels, density). keep=0 left/1 right."""
    if keep not in (0,1):
        raise ValueError('Invalid partial trace side')
    out = zero(16)
    for labels, matrix in blocks.values():
        for i, a in enumerate(labels):
            for j, b in enumerate(labels):
                if keep == 0 and a >> 4 == b >> 4:
                    out[a & 15][b & 15] += matrix[i][j]
                if keep == 1 and a & 15 == b & 15:
                    out[a >> 4][b >> 4] += matrix[i][j]
    return out


def verify_dual(payload):
    """Exact ceiling on lower bounds attainable by this boundary family.

    Compatible local densities need NOT extend to a physical global state.
    This output is not a physical ground-energy upper bound.
    """
    require_model(payload)
    raw = payload['dual_marginals']
    if not isinstance(raw, list) or len(raw) != 3:
        raise ValueError('Exactly three local marginals required')
    labels_by_key = sector_labels(4)
    keys = {f'{a},{b}' for a,b in labels_by_key}
    densities = []
    checked = 0
    for item in raw:
        if not isinstance(item, dict) or set(item) != keys:
            raise ValueError('All local charge blocks required exactly once')
        blocks = {}
        trace = F(0)
        for key, labels in labels_by_key.items():
            matrix = parse_matrix(item[f'{key[0]},{key[1]}'], len(labels))
            require_psd(matrix)
            trace += sum(matrix[i][i] for i in range(len(matrix)))
            blocks[key] = (labels, matrix)
            checked += 1
        if trace != 1:
            raise ValueError('Local marginal must have exact unit trace')
        densities.append(blocks)
    for i in range(2):
        if partial_trace(densities[i], 1) != partial_trace(densities[i+1], 0):
            raise ValueError('Adjacent marginal densities disagree on shared rung')
    bare = local_blocks([zero(16), zero(16)])
    ceiling = F(0)
    for rho, h in zip(densities, bare):
        for key, (labels, matrix) in rho.items():
            ham = h[key][1]
            ceiling += sum(matrix[i][j]*ham[j][i]
                           for i in range(len(labels)) for j in range(len(labels)))
    return {'status':'accepted_exact_family_ceiling',
            'family_lower_ceiling_over_t':str(ceiling),
            'local_sector_blocks_checked':checked,
            'global_determinants_enumerated':0,
            'physical_global_state_claimed':False,
            'scope':'Three adjacent-two-rung all-charge PSD blocks; arbitrary real symmetric spin-number-preserving one-rung boundary messages; fixed chemical potential 4'}
