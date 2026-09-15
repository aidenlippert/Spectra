"""Exact coefficient-space constructor for a diagonally dominant pair SOS.

This deliberately never builds a Fock-space matrix.  It parses a sparse
quartic pair Hamiltonian into a rational Gram dictionary, emits weighted
two-word squares, and verifies every coefficient exactly.
"""

from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import add, encode, scale, verify, word_product
from experiments.marginal_transfer_verify import apply_word

def car_check(modes):
    """Replay actual CAR H and integer square factors, not just a formal Gram."""
    labels, terms=family_hamiltonian(modes)
    word={u:((0,2*j+1),(0,2*j)) for j,u in enumerate(labels)}
    dagger=lambda w:tuple((1-c,i) for c,i in reversed(w))
    h=add(*(scale(dict(word_product(dagger(word[u]),word[v])),c) for u,v,c in terms))
    # Recognize H through the injective pair-word coordinate map. This explicit
    # implementation scans O(M^2) pairs; it makes no linear-time search claim.
    coordinates={}
    for u in labels:
        for v in labels:
            expansion=dict(word_product(dagger(word[u]),word[v]))
            assert len(expansion)==1
            w,sign=next(iter(expansion.items()))
            coordinates[w]=(u,v,sign)
    if set(h)-set(coordinates):raise ValueError('Outside the recognizable pair Hamiltonian family')
    parsed=[(*coordinates[w][:2],c/coordinates[w][2]) for w,c in h.items()]
    weighted=construct(parse_hamiltonian(labels,parsed))
    blocks=[]
    for j,(weight,vec) in enumerate(weighted):
        if weight.denominator!=1:raise ValueError('This integer demo requires integer weights')
        ws=[word[u] for u in vec]
        # Integer weights are repeated squares, avoiding an irrational sqrt.
        blocks.append({'name':f'pair-{j}','words':ws,'factor':[[int(v) for v in vec.values()]]*int(weight)})
    cert={'modes':modes,'particles':modes//2,'hamiltonian':encode(h),'b':'0',
          'number_multiplier':[],'denominator':1,'blocks':blocks}
    receipt=verify(cert)
    state=sum(1<<(2*j) for j in range(modes//2))
    assert all(apply_word(w,state) is None for w in word.values())
    assert receipt['lower']=='0' and receipt['residual_l1']=='0'
    return receipt


def construct(gram):
    labels = sorted(gram)
    cert = []
    for i, u in enumerate(labels):
        for v in labels[i + 1:]:
            x = gram[u].get(v, Fraction(0))
            if not x:
                continue
            cert.append((abs(x), {u: Fraction(1), v: Fraction(1 if x > 0 else -1)}))
    for u in labels:
        d = gram[u].get(u, Fraction(0)) - sum(
            abs(gram[u].get(v, Fraction(0))) for v in labels if v != u
        )
        if d < 0:
            raise ValueError(f"dominance fails for {u}: residual diagonal {d}")
        if d:
            cert.append((d, {u: Fraction(1)}))
    return cert


def expand(cert, labels):
    out = {u: {v: Fraction(0) for v in labels} for u in labels}
    for weight, vector in cert:
        for u, cu in vector.items():
            for v, cv in vector.items():
                out[u][v] += weight * cu * cv
    return out


def family_hamiltonian(modes):
    # Pair words are a_(2j+1) a_(2j), with sparse pair exchange on a chain.
    # This is a quartic pair-density/pair-hopping Hamiltonian family.
    labels = [f"p{j}" for j in range(modes // 2)]
    terms = []
    for j, u in enumerate(labels):
        terms.append((u, u, Fraction(3)))
        if j + 1 < len(labels):
            terms.append((u, labels[j + 1], Fraction(-1)))
            terms.append((labels[j + 1], u, Fraction(-1)))
    return labels, terms


def parse_hamiltonian(labels, terms):
    """Build G from H's sparse pair-word coefficients; reject asymmetry."""
    gram = {u: {} for u in labels}
    for u, v, coefficient in terms:
        if u not in gram or v not in gram:
            raise ValueError("unknown pair word")
        gram[u][v] = gram[u].get(v, Fraction(0)) + coefficient
    for u in labels:
        for v in labels:
            if gram[u].get(v, 0) != gram[v].get(u, 0):
                raise ValueError("Hamiltonian coefficient map is not Hermitian")
    return gram


def main():
    for modes in (8,16,32):
        receipt=car_check(modes)
        print(f'production CAR replay M={modes}: lower=upper=0, factor_rows={receipt["factor_rows"]}')
    labels, terms = family_hamiltonian(8)
    gram = parse_hamiltonian(labels, terms)
    cert = construct(gram)
    assert expand(cert, labels) == {u: {v: gram[u].get(v, 0) for v in labels} for u in labels}
    # Nontrivial exact lower bound: H = v^dagger G v is itself SOS, so b=0 is tight.
    assert len(cert) == 2 * len(labels) - 1
    # Half-filled witness |1,3,5,7> has at most one fermion in every pair,
    # so every pair-annihilation word kills it. Therefore H has exact energy 0.
    witness = {0, 2, 4, 6}
    assert all(not ({2 * j, 2 * j + 1} <= witness) for j in range(4))
    print(f"modes=8 pair_words={len(labels)} atoms={len(cert)} exact_lower_bound=0 tight_witness_N={len(witness)}")

    for modes in (16, 32):
        labels, terms = family_hamiltonian(modes)
        gram = parse_hamiltonian(labels, terms)
        cert = construct(gram)
        labels = sorted(gram)
        assert expand(cert, labels) == {u: {v: gram[u].get(v, 0) for v in labels} for u in labels}
        print(f"modes={modes} pair_words={len(labels)} atoms={len(cert)} global_fock_dim=2^{modes}")

    # PSD but non-dominant: direct recognizer rejects although a compact factor exists.
    bad = {"u": {"u": Fraction(1), "v": Fraction(2)},
           "v": {"u": Fraction(2), "v": Fraction(4)}}
    try:
        construct(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("non-dominant Gram must be rejected")
    print("PASS: exact constructive pair SOS, mode scaling, and rejection boundary")


if __name__ == "__main__":
    main()
