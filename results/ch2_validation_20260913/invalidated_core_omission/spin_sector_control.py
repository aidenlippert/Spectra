"""Small exact spin-sector control for the CH2 pilot.

This is a validation control, not a CH2 Hamiltonian.  It constructs the
two-electron/two-spatial-orbital determinant space and projects the M_S=0
block into S=0 and S=1 using S^2.  Keeping this separate makes it impossible
to accidentally label an M_S=0 calculation as a singlet calculation.
"""

from itertools import combinations
import numpy as np


def determinants(norb=2, nalpha=1, nbeta=1):
    alpha = list(combinations(range(norb), nalpha))
    beta = list(combinations(range(norb), nbeta))
    return [(a, b) for a in alpha for b in beta]


def spin_square_ms0(norb=2):
    """Return S^2 in the determinant basis, using a direct spin-flip action.

    For one alpha and one beta electron, S^2 = 1 + S_- S_+ in the M_S=0
    block.  The spin-flip term exchanges the two occupied spatial orbitals.
    """
    ds = determinants(norb)
    out = np.zeros((len(ds), len(ds)), dtype=float)
    index = {d: i for i, d in enumerate(ds)}
    for i, (a, b) in enumerate(ds):
        if a != b:
            out[i, i] = 1.0
            exchanged = ((b[0],), (a[0],))
            out[index[exchanged], i] += 1.0
    return ds, out


def exact_spin_gap():
    """A synthetic diagonal control whose S=0/S=1 gap is exactly 2 Ha.

    The diagonal operator assigns 0 to the singlet-coupled spatial pair and
    2 to the triplet-coupled pair after S^2 diagonalization.  It exists only
    to test state labeling and interval subtraction.
    """
    ds, s2 = spin_square_ms0()
    vals, vecs = np.linalg.eigh(s2)
    singlet = vecs[:, np.argmin(abs(vals - 0.0))]
    triplet = vecs[:, np.argmax(vals)]
    h = np.outer(triplet, triplet) * 2.0
    es = float(singlet @ h @ singlet)
    et = float(triplet @ h @ triplet)
    return {"singlet": es, "triplet": et, "gap_ET_minus_ES": et - es,
            "s2_eigenvalues": vals.tolist(), "determinants": ds}


if __name__ == "__main__":
    r = exact_spin_gap()
    assert all(abs(x-y) < 1e-12 for x, y in zip(r["s2_eigenvalues"], [0.0, 0.0, 0.0, 2.0]))
    assert abs(r["gap_ET_minus_ES"] - 2.0) < 1e-12
    print(r)
