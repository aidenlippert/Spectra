"""Reproducible numerical Lowdin orbital target generation."""
from fractions import Fraction as F
from itertools import product
from math import prod
import hashlib, json
import math
from pathlib import Path

from .marginal_symbolic import canonical, decode


def phase_alignment(generated, saved, spatial_orbitals):
    """Return the best phase pattern (first phase fixed positive) and errors."""
    if type(spatial_orbitals) is not int or not 1 <= spatial_orbitals <= 8:
        raise ValueError('spatial_orbitals must be an integer in 1..8')
    keys = set(generated) | set(saved)
    best = None
    for tail in product((-1, 1), repeat=spatial_orbitals - 1):
        signs = (1,) + tail
        diffs = [abs(float(generated.get(w, F(0))) *
                     prod(signs[i // 2] for _, i in w) -
                     float(saved.get(w, F(0)))) for w in keys]
        cand = (sum(diffs), max(diffs, default=0.0), signs)
        if best is None or cand[:2] < best[:2]:
            best = cand
    return best


def build_target(source, output, *, fixture_source=None, tolerance=1e-7):
    """Generate a numerical Lowdin target from fixture metadata and source H."""
    source = Path(source); output = Path(output)
    if output.exists(): raise ValueError('Refusing to overwrite existing target')
    data = json.loads(source.read_text())
    fixture = json.loads(Path(fixture_source).read_text()) if fixture_source else data
    modes = data.get('modes'); particles = data.get('particles')
    if type(modes) is not int or type(particles) is not int:
        raise ValueError('modes and particles must be integers')
    orbitals = modes // 2
    if modes < 2 or modes > 16 or modes % 2: raise ValueError('expected even modes 2..16')
    if particles < 0 or particles > modes: raise ValueError('particles out of range')
    geometry = fixture.get('geometry')
    if not geometry or not fixture.get('basis') or not fixture.get('unit'):
        raise ValueError('source lacks reproducible fixture geometry/basis/unit')
    if not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool) or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-6:
        raise ValueError('tolerance must be finite and in (0, 1e-6]')
    import numpy as np
    import pyscf
    from pyscf import gto, scf, ao2mo
    mol = gto.M(atom=geometry, basis=fixture['basis'], unit=fixture['unit'], spin=0,
                charge=0, verbose=0)
    if particles != mol.nelectron: raise ValueError('particle count does not match molecule')
    mf = scf.RHF(mol).run(conv_tol=1e-12, verbose=0)
    if not mf.converged or mol.nao_nr() != orbitals:
        raise ValueError('RHF convergence or orbital dimension mismatch')
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), orbitals).reshape((orbitals,) * 4)
    raw = {}
    for p in range(modes):
        for q in range(modes):
            v = float(h1[p//2, q//2]) if p % 2 == q % 2 else 0.
            if v: raw[((1,p),(0,q))] = F.from_float(v)
            for r in range(modes):
                for s in range(modes):
                    v = (eri[p//2,r//2,q//2,s//2]*(p%2==r%2)*(q%2==s%2)
                         - eri[p//2,s//2,q//2,r//2]*(p%2==s%2)*(q%2==r%2))/4
                    if v: raw[((1,p),(1,q),(0,s),(0,r))] = F.from_float(float(v))
    raw = canonical(raw)
    # The encoded fixture Hamiltonian has at most two-body (degree-four) words.
    saved = decode(data['hamiltonian'], modes, 4)
    align = phase_alignment(raw, saved, orbitals)
    if align[0] > tolerance: raise ValueError('canonical source phase alignment failed')
    overlap = mol.intor('int1e_ovlp'); vals, vec = np.linalg.eigh(overlap)
    half = (vec * np.sqrt(vals)) @ vec.T
    u = np.diag(align[2]) @ mf.mo_coeff.T @ half
    result = {'status':'aligned_numerical_proposal', 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'geometry':geometry, 'unit':fixture['unit'], 'basis':fixture['basis'], 'modes':modes, 'particles':particles,
      'pyscf_version':pyscf.__version__, 'orbital_basis':'Lowdin symmetric orthogonalized atomic orbitals',
      'source_orbital_basis':'RHF canonical MO',
      'rotation_convention':'U=diag(alignment_signs) C_MO.T S^(1/2); H_new=U.T H_old U; C_new=S^(-1/2)',
      'orthogonality_max_abs':float(np.max(np.abs(u.T@u-np.eye(orbitals)))),
      'canonical_phase_alignment':{'success':True,'phase_pattern':list(align[2]),'coefficient_l1_mismatch':align[0],
       'coefficient_max_abs_mismatch':align[1], 'generated_terms':len(raw),'saved_terms':len(saved), 'tolerance':tolerance},
      'target_rotation':u.tolist(),
      'scope':'Numerical aligned orbital target only; untrusted proposal, with exact rational replacement and Hamiltonian rotation verified separately.'}
    output.write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', required=True); ap.add_argument('--fixture-source')
    ap.add_argument('--output', required=True); ap.add_argument('--tolerance', type=float, default=1e-7)
    args = ap.parse_args()
    print(json.dumps(build_target(args.source, args.output, fixture_source=args.fixture_source,
                                  tolerance=args.tolerance), indent=2))
