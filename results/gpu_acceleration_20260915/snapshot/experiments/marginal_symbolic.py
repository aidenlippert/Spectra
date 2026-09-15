"""CAR certificate extraction and exact symbolic lower-bound verification.

The verifier uses sparse polynomials, not occupation-sector matrices. NumPy is
loaded only by the optional least-squares proposer. All accepted bounds come
from exact rational re-expansion and a coefficient l1 residual bound.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
import argparse
import hashlib
import json
from pathlib import Path

from experiments.marginal_hunt_car import add, adj, mono, mul, scale


@lru_cache(maxsize=100000)
def word_product(left, right):
    return tuple(mul(mono(left), mono(right)).items())


def product(p, q):
    result = {}
    for left, a in p.items():
        for right, b in q.items():
            for word, sign in word_product(left, right):
                result[word] = result.get(word, 0) + a*b*sign
    return {w: c for w, c in result.items() if c}


def canonical(p):
    # Boolean creation flags historically came from adj(). They compare and
    # hash equal to integers, so cached products can otherwise leak boolean
    # labels into a later strictly validated integer-word calculation.
    return {tuple((int(c) if type(c) is bool else c,i) for c,i in w):v
            for w,v in product(p, {(): 1}).items()}


def hermitian(p):
    return canonical(adj(p)) == canonical(p)


def number_shift(modes, particles):
    return add(mono((), -particles), *(mono(((1, i), (0, i))) for i in range(modes)))


def validate_word(word, modes, degree):
    if len(word) > degree or any(len(letter) != 2 for letter in word):
        raise ValueError("Invalid word length")
    if any(type(c) is not int or c not in (0, 1) or type(i) is not int or not 0 <= i < modes for c, i in word):
        raise ValueError("Invalid ladder operator")
    return tuple(tuple(letter) for letter in word)


def expand_squares(blocks, denominator, modes, max_degree=3):
    if type(denominator) is not int or denominator <= 0:
        raise ValueError("Invalid factor denominator")
    total = {}
    rows = 0
    nonzeros = 0
    for block in blocks:
        words = [validate_word(w, modes, max_degree) for w in block["words"]]
        charges = {sum(2*c-1 for c, _ in w) for w in words}
        if not words or len(charges) != 1:
            raise ValueError("Each factor dictionary must have one charge")
        factors = block["factor"]
        if any(len(row) != len(words) or any(type(c) is not int for c in row) for row in factors):
            raise ValueError("Invalid integer factor")
        rows += len(factors)
        nonzeros += sum(c != 0 for row in factors for c in row)
        # Exact Gram coefficients avoid normal-ordering each dense square anew.
        # The real Gram matrix is symmetric. Compute each dot product once,
        # with contiguous columns, and emit both ordered CAR products.
        columns = [tuple(row[i] for row in factors) for i in range(len(words))]
        adjoints = [tuple((1-c, m) for c, m in reversed(w)) for w in words]
        for i, wi in enumerate(words):
            for j in range(i, len(words)):
                coefficient = sum(a*b for a,b in zip(columns[i], columns[j]))
                if not coefficient:
                    continue
                for word, sign in word_product(adjoints[i], words[j]):
                    total[word] = total.get(word, 0) + coefficient*sign
                if i != j:
                    for word, sign in word_product(adjoints[j], words[i]):
                        total[word] = total.get(word, 0) + coefficient*sign
    total = {w: F(c, denominator**2) for w, c in total.items() if c}
    return total, {"factor_rows": rows, "factor_nonzeros": nonzeros}


def encode(poly):
    return [{"word": [[int(c), i] for c, i in w], "coefficient": str(c)}
            for w, c in sorted(poly.items(), key=lambda x: (len(x[0]), x[0])) if c]


def decode(terms, modes, max_degree):
    raw = {}
    for term in terms:
        word = validate_word(term["word"], modes, max_degree)
        value = term["coefficient"]
        if not isinstance(value, str):
            raise ValueError("Coefficients must be exact rational strings")
        coefficient = F(value)
        raw[word] = raw.get(word, 0) + coefficient
    return canonical(raw)


def transform(poly, mapping):
    return canonical({tuple((c, mapping[i]) for c, i in w): a for w, a in poly.items()})


def signed_orbit(word, mappings):
    orbit = {}
    for mapping in mappings:
        for image, sign in transform(mono(word), mapping).items():
            if image in orbit and orbit[image] != sign:
                raise ValueError("Orbit contains inconsistent stabilizer signs")
            orbit[image] = sign
    return orbit


def expand_orbits(recipe, modes):
    mappings = recipe["permutations"]
    if not mappings or any(any(type(i) is not int for i in p) or sorted(p) != list(range(modes)) for p in mappings):
        raise ValueError("Invalid orbital permutation")
    terms = []
    for item in recipe["orbits"]:
        word = validate_word(item["representative"], modes, 4)
        if not isinstance(item["coefficient"], str):
            raise ValueError("Orbit coefficient must be a rational string")
        terms.append(scale(signed_orbit(word, mappings), F(item["coefficient"])))
    return add(*terms)


def verified_residual(certificate):
    """Return the validated exact residual and its coefficient-norm receipt.

    Residual-bound extensions reuse the same validation and CAR expansion;
    they must not bypass this function's input and Hermiticity gates.
    """
    modes, particles = certificate["modes"], certificate["particles"]
    if type(modes) is not int or modes < 1 or type(particles) is not int or not 0 <= particles <= modes:
        raise ValueError("Invalid fermion sector")
    degree = certificate.get("operator_degree", 3)
    if type(degree) is not int or degree not in (3, 4):
        raise ValueError("Supported factor degrees are 3 and 4")
    h = decode(certificate["hamiltonian"], modes, 4)
    if ("number_multiplier" in certificate) == ("multiplier_orbits" in certificate):
        raise ValueError("Provide exactly one multiplier representation")
    x = (decode(certificate["number_multiplier"], modes, 2*degree-2) if "number_multiplier" in certificate
         else expand_orbits(certificate["multiplier_orbits"], modes))
    for p in (h, x):
        if any(sum(2*c-1 for c, _ in w) != 0 for w in p) or not hermitian(p):
            raise ValueError("H and X must be Hermitian and number conserving")
    if not isinstance(certificate["b"], str):
        raise ValueError("b must be a rational string")
    b = F(certificate["b"])
    squares, stats = expand_squares(certificate["blocks"], certificate["denominator"], modes, degree)
    ideal = product(number_shift(modes, particles), x)
    residual = add(h, mono((), -b), scale(squares, -1), scale(ideal, -1))
    if not hermitian(residual):
        raise ValueError("Non-Hermitian residual")
    eta = sum(abs(c) for c in residual.values())
    # Every normal ordered monomial is a product of contractions, hence norm<=1.
    receipt = {"lower": str(b-eta), "lower_float": float(b-eta),
            "residual_l1": str(eta), "residual_l1_float": float(eta),
            "multiplier_terms": len(x), "square_polynomial_terms": len(squares),
            "residual_terms": len(residual), "residual_max_degree": max(map(len, residual), default=0),
            "mode_count": modes, "particle_number": particles, **stats}
    return residual, receipt


def verify(certificate):
    """Sound for any admitted finite M,N; no sector enumeration or eigenvalues."""
    return verified_residual(certificate)[1]


def multiplier_basis(modes, max_body=2):
    basis = []
    for k in range(max_body+1):
        sets = list(combinations(range(modes), k))
        for i, left in enumerate(sets):
            for right in sets[i:]:
                word = tuple((1, p) for p in left) + tuple((0, p) for p in reversed(right))
                e = canonical(mono(word))
                basis.append(e if left == right else scale(add(e, canonical(adj(e))), F(1, 2)))
    return basis


def model_hamiltonian(t, asymmetric=False):
    """Polynomial form of the prior six-mode hopping test, built without Fock matrices."""
    t = F(t)
    n = [mono(((1, i), (0, i))) for i in range(6)]
    h = add(*(product(n[i], n[j]) for triple in ((0, 1, 2), (3, 4, 5)) for i, j in combinations(triple, 2)))
    edges = [(0, 3, F(1)), (1, 4, F(1)), (2, 5, F(1))]
    if asymmetric:
        edges = [(0, 3, F(1)), (1, 4, F(7, 10)), (2, 5, F(13, 10)), (0, 4, F(2, 5))]
    for i, j, weight in edges:
        e = mono(((1, i), (0, j)))
        h = add(h, scale(add(e, canonical(adj(e))), -t*weight))
    return h


def extract(original, rounding=10**10):
    """Numerically proposes X; verify() alone decides the sound lower bound."""
    import numpy as np
    modes, particles = original["mode_count"], original["particle_number"]
    if modes != 6 or particles != 3:
        raise ValueError("This importer targets the prior six-mode artifacts")
    h = model_hamiltonian(original["t"], original["asymmetric"])
    b = F(original["b_numerator"], original["denominator"]**2)
    squares, _ = expand_squares(original["blocks"], original["denominator"], modes)
    target = add(h, mono((), -b), scale(squares, -1))
    basis = multiplier_basis(modes)
    shift = number_shift(modes, particles)
    columns = [product(shift, q) for q in basis]
    words = sorted(set(target).union(*(set(q) for q in columns)), key=lambda w: (len(w), w))
    matrix = np.array([[float(col.get(w, 0)) for col in columns] for w in words])
    rhs = np.array([float(target.get(w, 0)) for w in words])
    coefficients, _, rank, _ = np.linalg.lstsq(matrix, rhs, rcond=None)
    rational = [F(int(round(c*rounding)), rounding) for c in coefficients]
    x = add(*(scale(q, c) for q, c in zip(basis, rational)))
    certificate = {"modes": modes, "particles": particles, "hamiltonian": encode(h),
                   "number_multiplier": encode(x), "b": str(b),
                   "denominator": original["denominator"], "blocks": original["blocks"]}
    receipt = verify(certificate)
    receipt.update({"raw_residual_l1": str(sum(abs(c) for c in target.values())),
                    "multiplier_basis_dimension": len(basis), "numeric_map_rank": int(rank),
                    "coefficient_rows": len(words), "multiplier_rounding_denominator": rounding})
    return certificate, receipt


def compress_matched(certificate):
    """Six-mode symmetry proposal; exact verification charges all compression error."""
    if certificate["modes"] != 6 or "number_multiplier" not in certificate:
        raise ValueError("Expected an uncompressed six-mode certificate")
    x = decode(certificate["number_multiplier"], 6, 4)
    mappings = [tuple(p[i % 3] + 3*((i // 3) ^ swap) for i in range(6))
                for p in permutations(range(3)) for swap in (0, 1)]
    sym = scale(add(*(transform(x, p) for p in mappings)), F(1, len(mappings)))
    sym = {w: a for w, a in sym.items() if abs(a) >= F(1, 1000000)}
    unseen, orbits = set(sym), []
    while unseen:
        w = min(unseen, key=lambda w: (len(w), w))
        orbit = signed_orbit(w, mappings)
        orbits.append({"representative": [[int(c), i] for c, i in w],
                       "coefficient": str(sym[w]), "orbit_size": len(orbit)})
        unseen -= set(orbit)
    recipe = {"permutations": [list(p) for p in mappings], "orbits": orbits}
    if expand_orbits(recipe, 6) != sym:
        raise AssertionError("Orbit recipe does not reproduce the multiplier")
    result = dict(certificate)
    del result["number_multiplier"]
    result["multiplier_orbits"] = recipe
    return result, verify(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--compress", type=Path)
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(verify(json.loads(args.verify.read_text())), indent=2))
        return
    if args.compress:
        certificate, receipt = compress_matched(json.loads(args.compress.read_text()))
        destination = args.compress.parent / "matched_symmetry_compressed.json"
        destination.write_text(json.dumps(certificate) + "\n")
        (args.compress.parent / "matched_symmetry_orbits.json").write_text(json.dumps(
            {**certificate["multiplier_orbits"], "receipt": receipt}, indent=2) + "\n")
        print(json.dumps(receipt, indent=2))
        return
    root = Path(__file__).resolve().parents[1]
    sources = [args.source] if args.source else [root / "results/marginal_hopping" / name for name in
        ("full_mixed_matched.json", "full_mixed_asymmetric.json", "matched_1_paired_mixed.json", "asymmetric_1_paired_mixed.json")]
    out = root / "results/marginal_symbolic"
    out.mkdir(exist_ok=True)
    summary = []
    for path in sources:
        source_bytes = path.read_bytes()
        certificate, receipt = extract(json.loads(source_bytes))
        certificate["source_sha256"] = hashlib.sha256(source_bytes).hexdigest()
        destination = out / path.name
        destination.write_text(json.dumps(certificate) + "\n")
        receipt.update({"source": str(path), "certificate_file": destination.name,
                        "source_sha256": certificate["source_sha256"]})
        summary.append(receipt)
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
