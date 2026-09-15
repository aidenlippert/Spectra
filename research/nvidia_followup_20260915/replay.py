"""Complete original-checker replay, including the H10 orbital-rotation upper."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
from research.transfer_solver_20260915.actions import run as original_actions
from research.transfer_solver_20260915.budget import dump
from research.correlated_pair_20260913.mps_exact import check as check_MPS
from research.molecular_collective_20260913.core import digest


def enable_flint_rationals():
    import flint
    from fractions import Fraction as PythonFraction
    # Load the complete lower-checker dependency graph before choosing its
    # arithmetic, and discard coefficient objects cached by the upper check.
    from research.collective_completion_20260914 import spin_screen
    from experiments import marginal_symbolic
    for name,module in list(sys.modules.items()):
        if name.startswith(('research.','experiments.')):
            for alias in ('F','Fraction'):
                if getattr(module,alias,None) is PythonFraction:
                    setattr(module,alias,flint.fmpq)
    marginal_symbolic.word_product.cache_clear()


def upper(case, rotated=None):
    read=lambda path:json.loads(path.read_text())
    data=read(case/'fixture.json')
    supplied=read(case/'upper.json')
    if rotated is None:
        checked=check_MPS(data,read(case/'mps/state.json'))
        value=F(checked['upper_Ha'])
    else:
        from research.transfer_followup_20260915.rotated_upper import rotate_and_round
        rotation=read(rotated/'rotation.json')
        if rotation['original_fixture_sha256']!=digest(data):
            raise ValueError('Rotation belongs to a different original Hamiltonian')
        recomputed, allowance=rotate_and_round(data,rotation['integer_matrix'],int(rotation['denominator']))
        if recomputed!=read(rotated/'fixture.json'):
            raise ValueError('Exact transformed coefficients do not reproduce')
        state=read(rotated/'mps/state.json')
        checked=check_MPS(recomputed,state)
        value=F(checked['upper_Ha'])+allowance
        checked={'upper_Ha':str(value),'rotated_MPS_upper':checked,
            'rotation_rounding_allowance_Ha':str(allowance),'original_fixture_sha256':digest(data),
            'exact_orthogonality_rechecked':True}
    if value!=F(supplied['upper_Ha']):
        raise ValueError('Upper endpoint does not reproduce from its actual witness')
    return checked


def run(case, proposal, output, rotated=None, compiled_grams=False, compiled_rationals=False):
    start=time.monotonic()
    output.mkdir(parents=True,exist_ok=False)
    receipt=upper(case,rotated)
    upper_seconds=time.monotonic()-start
    for name in ('fixture.json','nonsinglet.json'):
        shutil.copyfile(case/name,output/name)
    dump(output/'upper.json',receipt)
    source=case/proposal/'export/certificate.json'
    target=output/'stage2/export/certificate.json'
    target.parent.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(source,target)
    if compiled_grams:
        from experiments import marginal_symbolic
        from research.nvidia_followup_20260915.flint_squares import expand_squares
        marginal_symbolic.expand_squares=expand_squares
    if compiled_rationals:
        enable_flint_rationals()
    original_actions(output,'replay')
    forbidden=[n for n in ('numpy','scipy','numba','cupy','quimb','pyscf','cvxpy') if n in sys.modules]
    if forbidden:
        raise ValueError(('Numerical accepting import',forbidden))
    dump(output/'complete_replay.json',{'proposal':str(source),
        'proposal_file_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'upper_seconds':upper_seconds,'complete_seconds':time.monotonic()-start,
        'all_upper_and_lower_dependencies_rechecked':True,
        'integer_Gram_backend':'FLINT' if compiled_grams else 'Python_standard_library',
        'rational_backend':'FLINT' if compiled_rationals else 'Python_standard_library',
        'independent_stdlib_comparison_required_for_experimental_backend':compiled_grams or compiled_rationals,
        'rotated_state_directory':str(rotated) if rotated else None,'numerical_imports':forbidden})


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('case',type=Path)
    p.add_argument('proposal')
    p.add_argument('output',type=Path)
    p.add_argument('--rotated',type=Path)
    p.add_argument('--compiled-grams',action='store_true')
    p.add_argument('--compiled-rationals',action='store_true')
    a=p.parse_args()
    run(a.case.resolve(),a.proposal,a.output.resolve(),a.rotated.resolve() if a.rotated else None,a.compiled_grams,a.compiled_rationals)
