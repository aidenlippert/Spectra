"""Exact compiled integer Gram products with the original CAR expansion gates."""
from fractions import Fraction as F
from experiments.marginal_symbolic import validate_word, word_product


def expand_squares(blocks, denominator, modes, max_degree=3):
    from flint import fmpz_mat
    if type(denominator) is not int or denominator <= 0:
        raise ValueError('Invalid factor denominator')
    total = {}
    rows = nonzeros = 0
    for block in blocks:
        words = [validate_word(w, modes, max_degree) for w in block['words']]
        charges = {sum(2*c-1 for c, _ in w) for w in words}
        if not words or len(charges) != 1:
            raise ValueError('Each factor dictionary must have one charge')
        factors = block['factor']
        if any(len(row) != len(words) or any(type(c) is not int for c in row) for row in factors):
            raise ValueError('Invalid integer factor')
        rows += len(factors)
        nonzeros += sum(c != 0 for row in factors for c in row)
        matrix = fmpz_mat(factors) if factors else fmpz_mat(0, len(words))
        gram = matrix.transpose()*matrix
        adjoints = [tuple((1-c, m) for c, m in reversed(w)) for w in words]
        for i in range(len(words)):
            for j in range(i, len(words)):
                coefficient = int(gram[i, j])
                if not coefficient:
                    continue
                for word, sign in word_product(adjoints[i], words[j]):
                    total[word] = total.get(word, 0)+coefficient*sign
                if i != j:
                    for word, sign in word_product(adjoints[j], words[i]):
                        total[word] = total.get(word, 0)+coefficient*sign
    return ({w: F(c, denominator**2) for w, c in total.items() if c},
            {'factor_rows': rows, 'factor_nonzeros': nonzeros})
