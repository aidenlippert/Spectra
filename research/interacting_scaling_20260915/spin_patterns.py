"""Reuse exact SU(2) coefficient patterns under ordered spatial relabeling."""
from fractions import Fraction as F
from functools import lru_cache
from research.sector_quotient_20260914.fast_twirl import twirl as original_twirl


@lru_cache(maxsize=16384)
def pattern(word):
    return tuple(original_twirl({word: F(1)}).items())


def twirl(poly):
    out = {}
    for word, value in poly.items():
        creators = tuple((c, i) for c, i in word if c)
        annihilators = tuple((c, i) for c, i in word if not c)
        if (len(word) > 6 or len(word) % 2 or word != creators+annihilators or
            creators != tuple(sorted(set(creators))) or annihilators != tuple(sorted(set(annihilators))) or
            any(c not in (0, 1) or type(i) is not int or i < 0 for c, i in word) or type(value) not in (int, F)):
            raise ValueError('Exact even canonical CAR polynomial through degree six required')
        spatial = sorted({i//2 for c, i in word})
        index = {i: j for j, i in enumerate(spatial)}
        key = tuple((c, 2*index[i//2]+i%2) for c, i in word)
        for local_word, coefficient in pattern(key):
            target = tuple((c, 2*spatial[i//2]+i%2) for c, i in local_word)
            out[target] = out.get(target, F(0))+value*coefficient
    return {w: value for w, value in out.items() if value}
