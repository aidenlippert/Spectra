"""Sparse determinant actions and exact branch proofs of complement row bounds."""
from fractions import Fraction as F
from math import comb

from experiments.marginal_symbolic import decode, hermitian
from experiments.marginal_transfer_verify import apply_word


class GapFailure(ValueError):
    def __init__(self, state, lower):
        super().__init__('A complementary determinant fails the requested row bound')
        self.state, self.lower = state, lower


class DeterminantOracle:
    def __init__(self, certificate):
        self.modes, self.particles = certificate.get('modes'), certificate.get('particles')
        if (type(self.modes) is not int or not 2 <= self.modes <= 64
                or type(self.particles) is not int or not 0 < self.particles < self.modes):
            raise ValueError('Expected 2..64 modes and an interior fixed particle number')
        self.h = decode(certificate['hamiltonian'], self.modes, 4)
        if not hermitian(self.h) or any(sum(2 * c - 1 for c, _ in word) for word in self.h):
            raise ValueError('Real Hermitian number-conserving Hamiltonian required')
        self.all_bits = (1 << self.modes) - 1
        self.cache, self.diagonal, self.transitions = {}, {}, []
        for word, coefficient in self.h.items():
            required, occupied, flip, parity = 0, 0, 0, 0
            for creation, mode in reversed(word):
                bit = 1 << mode
                need = (1 - creation) ^ bool(flip & bit)
                if required & bit and bool(occupied & bit) != need:
                    break
                required |= bit
                if need:
                    occupied |= bit
                if (flip & (bit - 1)).bit_count() % 2:
                    coefficient = -coefficient
                parity ^= bit - 1
                flip ^= bit
            else:
                if flip:
                    self.transitions.append((required, occupied, flip, abs(coefficient)))
                else:
                    if required != occupied or parity & ~required or required.bit_count() > 2:
                        raise ValueError('Diagonal must reduce to a quadratic occupation polynomial')
                    value = coefficient * (-1 if (occupied & parity).bit_count() % 2 else 1)
                    self.diagonal[required] = self.diagonal.get(required, F(0)) + value

    def valid_state(self, state):
        return type(state) is int and 0 <= state <= self.all_bits and state.bit_count() == self.particles

    @property
    def sector_dimension(self):
        return comb(self.modes, self.particles)

    def retained(self, states):
        if (type(states) is not list or not 1 <= len(states) <= 32
                or len(set(states)) != len(states) or any(not self.valid_state(s) for s in states)):
            raise ValueError('Expected one to 32 distinct physical retained determinants')
        if len(states) == self.sector_dimension:
            raise ValueError('A nonempty complementary sector is required')
        return set(states)

    def action(self, state):
        if not self.valid_state(state):
            raise ValueError('Invalid determinant')
        if state not in self.cache:
            if len(self.cache) >= 4096:
                raise ValueError('Sparse action budget exhausted')
            result = {}
            for word, coefficient in self.h.items():
                image = apply_word(word, state)
                if image:
                    target, sign = image
                    result[target] = result.get(target, F(0)) + sign * coefficient
            self.cache[state] = {s: value for s, value in result.items() if value}
        return self.cache[state]

    def referenced_state_count(self):
        return len(set(self.cache) | {s for row in self.cache.values() for s in row})

    def diagonal_lower(self, mask, bits):
        """Fixed-N lower: k smallest linear terms and choose(k,2) pair terms."""
        base = F(0)
        linear = {i: F(0) for i in range(self.modes) if not mask & (1 << i)}
        pairs = {}
        for support, coefficient in self.diagonal.items():
            if support & mask & ~bits:
                continue
            free = support & ~mask
            if not free:
                base += coefficient
            elif free.bit_count() == 1:
                linear[free.bit_length() - 1] += coefficient
            else:
                pairs[free] = pairs.get(free, F(0)) + coefficient
        k = self.particles - bits.bit_count()
        pair_values = [pairs.get((1 << i) | (1 << j), F(0)) for i in linear for j in linear if i < j]
        return base + sum(sorted(linear.values())[:k]) + sum(sorted(pair_values)[:k * (k - 1) // 2])

    def branch_lower(self, mask, bits, retained):
        lower = self.diagonal_lower(mask, bits)
        for required, occupied, flip, magnitude in self.transitions:
            if (bits ^ occupied) & mask & required:
                continue
            fixed, source = mask | required, bits | occupied
            free = self.modes - fixed.bit_count()
            k = self.particles - source.bit_count()
            if not 0 <= k <= free:
                continue
            possible = comb(free, k)
            # Count the inverse images of P, not completions of the Q branch.
            excluded = sum(1 for p in retained if ((p ^ flip) & fixed) == source
                           and (p ^ flip).bit_count() == self.particles)
            if excluded != possible:
                lower -= magnitude
        return lower

    def cover(self, states, gamma, tree=None):
        retained = self.retained(states)
        gamma = F(gamma)
        building = tree is None
        stats = {'nodes': 0, 'pruned_branches': 0, 'exact_Q_rows': 0,
                 'pruned_Q_states': 0, 'covered_sector_states': 0}

        def visit(mask, bits, node):
            stats['nodes'] += 1
            if stats['nodes'] > 100000:
                raise ValueError('Complement proof node budget exhausted')
            free = self.modes - mask.bit_count()
            k = self.particles - bits.bit_count()
            if not 0 <= k <= free:
                raise ValueError('Invalid particle-count branch')
            if k == 0:
                mask = self.all_bits
            elif k == free:
                bits |= self.all_bits ^ mask
                mask = self.all_bits
            if not building and (type(node) is not list or not node or type(node[0]) is not str):
                raise ValueError('Malformed complement proof node')
            if mask == self.all_bits:
                stats['covered_sector_states'] += 1
                expected = ['retained'] if bits in retained else ['row']
                if not building and node != expected:
                    raise ValueError('Incorrect terminal complement proof node')
                if bits not in retained:
                    stats['exact_Q_rows'] += 1
                    row = self.action(bits)
                    lower = row.get(bits, F(0)) - sum(abs(x) for s, x in row.items() if s != bits and s not in retained)
                    if lower < gamma:
                        raise GapFailure(bits, lower)
                return expected
            lower = self.branch_lower(mask, bits, retained)
            if (building and lower >= gamma) or (not building and node == ['bound']):
                if lower < gamma:
                    raise ValueError('Unproved partial-occupation bound')
                count = comb(free, k)
                stats['pruned_branches'] += 1
                stats['covered_sector_states'] += count
                stats['pruned_Q_states'] += count - sum(1 for p in retained if p & mask == bits)
                return ['bound']
            if not building and (len(node) != 3 or node[0] != 'split'):
                raise ValueError('A split must include both complementary branches')
            bit = 1 << next(i for i in range(self.modes) if not mask & (1 << i))
            left = visit(mask | bit, bits, None if building else node[1])
            right = visit(mask | bit, bits | bit, None if building else node[2])
            return ['split', left, right]

        result = visit(0, 0, tree)
        if stats['covered_sector_states'] != comb(self.modes, self.particles):
            raise ValueError('Complement proof failed complete fixed-number coverage')
        return result, stats

    def retained_data(self, states, gamma):
        retained = self.retained(states)
        columns = [self.action(s) for s in states]
        block = [[columns[j].get(s, F(0)) for j in range(len(states))] for s in states]
        leakage = [[sum(value * columns[j].get(s, F(0)) for s, value in column.items() if s not in retained)
                    for j in range(len(states))] for column in columns]
        return block, leakage, F(gamma)

    def upper(self, witness):
        if type(witness) is not dict:
            raise ValueError('Sparse integer upper witness required')
        states, amplitudes = witness.get('states'), witness.get('amplitudes')
        if (type(states) is not list or not 1 <= len(states) <= 4096
                or len(set(states)) != len(states) or any(not self.valid_state(s) for s in states)
                or type(amplitudes) is not list or len(states) != len(amplitudes)
                or any(type(x) is not int for x in amplitudes) or not any(amplitudes)):
            raise ValueError('Invalid sparse integer upper witness')
        vector = dict(zip(states, amplitudes))
        norm = sum(x * x for x in amplitudes)
        energy = sum(x * y * self.action(s).get(t, F(0)) for s, x in vector.items() if x for t, y in vector.items() if y)
        return energy / norm
