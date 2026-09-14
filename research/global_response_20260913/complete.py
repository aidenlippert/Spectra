"""Integrated exact replay with explicit inherited terminal dependencies.

This completes reference-assisted certificates, not new cheap SOS discovery.
The direct fresh-case procedure is replayed even when its interval is wide.
Default replay writes nothing; --build creates this campaign's derived files.
"""
from copy import deepcopy
from fractions import Fraction as F
import argparse
import hashlib
import json
import resource
import sys
import time

from experiments.marginal_symbolic import decode, verified_residual
from research.certificate_scaling.streaming_reference_upper import upper
from research.certificate_scaling import wedge_spectral_bound as spectral
from research.compact_response_20260913 import program, closure
from research.composable_response_20260913 import joint
from research.ch2_validation_20260913 import certify as ch2
from research.global_response_20260913 import global_program as g, slater
from research.global_response_20260913.reference_diagnostic import lower, scalar_terminal_margin
from research.global_response_20260913.lifted_upper import lift


def read(path):
    return json.loads(path.read_text())


def file_binding(path):
    raw = path.read_bytes()
    return {'path': str(path.relative_to(g.ROOT)), 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}


def bind_hamiltonian(data, cert):
    if (data['modes'], data['particles']) != (cert['modes'], cert['particles']) or decode(data['hamiltonian'], data['modes'], 4) != decode(cert['hamiltonian'], cert['modes'], 4):
        raise ValueError('Hamiltonian or physical-sector mismatch')


def response_variant(data, tail, original, target, budget, path, build):
    if build:
        cert = deepcopy(original)
        cert['target_Ha'], cert['budget_Ha'] = str(target), str(budget)
        bound = g.bounds(data, tail, target, cert['first_sector'], cert['joint_sector'])
        cert['partition'] = bound['partition']; cert['response'] = g.scalar_program(bound, budget)
        path.write_text(json.dumps(cert, indent=2)+'\n')
    cert = read(path)
    if F(cert['target_Ha']) != target or F(cert['budget_Ha']) != budget:
        raise ValueError('Terminal-aware recipe changed')
    return cert, g.check(data, tail, cert)


def reference_closed(name, build):
    data, tail, _ = g.load_case(name)
    L, inherited = lower(data, name)
    meta = read(g.ROOT/f'results/certificate_scaling/cubic_precision/intervals/final_{name}.json')
    ref = read(g.ROOT/meta['reference']); U, ustats = upper(data, ref['independent_upper'])
    if U != F(meta['upper']): raise ValueError('Reference upper did not reproduce')
    original = read(g.OUT/f'{name}_global.json'); response = g.check(data, tail, original)
    b = F(original['target_Ha']); gamma = L-b
    if gamma <= 0: raise ValueError('Inherited terminal margin unavailable')
    margin = scalar_terminal_margin(L, b, [F(response['residual_penalty_Ha'])])
    if not margin['positive']: raise ValueError('Uniform response does not close')
    budget = min(gamma/2, F(1, 1000))
    cert, tuned = response_variant(data, tail, original, b, budget, g.OUT/f'{name}_terminal_tuned.json', build)
    tuned_margin = scalar_terminal_margin(L, b, [F(tuned['residual_penalty_Ha'])])
    if not tuned_margin['positive'] or U < b: raise ValueError('Integrated positivity failed')
    result = {'case': name, 'lower_Ha': str(b), 'upper_Ha': str(U), 'width_mHa': float((U-b)*1000),
              'stronger_inherited_lower_Ha': str(L), 'uniform_response': response,
              'uniform_terminal': margin, 'terminal_aware_response': tuned,
              'terminal_aware_margin': tuned_margin, 'inherited_terminal_proof': inherited,
              'upper_replay': ustats, 'upper_dependency': file_binding(g.ROOT/meta['reference']),
              'new_small_terminal_discovery': False, 'full_fixed_N_enumeration_in_lower_replay': False,
              'scope': 'The inherited cubic SOS proof supplies terminal positivity. Response tuning spends at most half its verified margin; it does not improve the inherited lower.'}
    if name == 'h8':
        _, _, old_b = g.load_case('h8_legacy_target')
        if U >= old_b: raise ValueError('Old-target counterexample did not reproduce')
        result['old_target_negative_expectation_Ha'] = str(U-old_b)
    return result


def transferred_closure(build):
    data, tail, _ = g.load_case('fresh_h6_1p6')
    old = read(joint.OUT/'transfer_program.json'); proof = read(joint.OUT/'transfer_closure.json')
    ref = read(joint.OUT/'transfer_upper.json')
    receipt = closure.check(data, tail, old, proof, ref)
    L, U = F(receipt['lower_Ha']), F(receipt['upper_Ha'])
    # A controlled 0.010 mHa loss supplies an explicit global margin for
    # transporting the inherited closure to this different response.
    b = L-F(1, 100000); original = read(g.OUT/'fresh_h6_1p6_global.json')
    _, response = response_variant(data, tail, original, b, F(1, 10**6), g.OUT/'fresh_h6_1p6_closed_global.json', build)
    terminal = scalar_terminal_margin(L, b, [F(response['residual_penalty_Ha'])])
    if not terminal['positive']: raise ValueError('Transferred inherited closure did not transport')
    return {'case': 'fresh_h6_1p6', 'preserved_width_mHa': receipt['width_mHa'],
            'new_reference_assisted_width_mHa': float((U-b)*1000), 'lower_Ha': str(b), 'upper_Ha': str(U),
            'response': response, 'terminal': terminal, 'inherited_closure': receipt,
            'scope': 'The transported closure still enumerates 924 fixed-N labels. No compact terminal replacement is claimed.'}


def fresh_full_interval():
    data, tail, _ = g.load_case('fresh_h6_1p73')
    cert = read(g.OUT/'fresh_h6_1p73_global.json'); response = g.check(data, tail, cert)
    path = g.OUT/'joint_factor/fresh_h6_1p73/certificate.json'; direct = read(path)
    bind_hamiltonian(data, direct)
    residual, base = verified_residual(direct)
    if max(map(len, residual), default=0) > 6: raise ValueError('Higher-body residual is outside this extension')
    # Candidate factors can come from numerical analysis, but acceptance is
    # against the ORIGINAL unmodified certificate and independently expanded
    # residual. No changed operator-degree tag is admitted here.
    proofpath = g.OUT/'joint_factor/spectral_local/fresh_h6_1p73/witness.json'
    sr = spectral.replay_residual(direct, residual, read(proofpath))
    ref = read(g.OUT/'fresh_h6_1p73/reference_upper.json'); U, us = upper(data, ref['independent_upper'])
    L = F(sr['lower'])
    if U < L: raise ValueError('Fresh upper/lower conflict')
    return {'case': 'fresh_h6_1p73', 'response': response, 'lower_Ha': str(L), 'upper_Ha': str(U),
            'width_Ha': str(U-L), 'width_mHa': float((U-L)*1000), 'passes_1p6mHa': U-L <= F(16, 10000),
            'direct_original_certificate': file_binding(path), 'raw_b_Ha': direct['b'], 'coefficient_replay': base,
            'spectral_replay': sr, 'candidate_factor_source': file_binding(proofpath),
            'accepted_against_unchanged_original_certificate': True, 'reference_upper_replay': us,
            'scope': 'Full interval for the frozen direct lower procedure. The response component passes, but this lower does not close its tight terminal target.'}


def compact_uppers():
    rows = []
    for row in read(g.OUT/'upper_model/compact_certificates.json'):
        data = read(g.ROOT/row['fixture']); r = slater.check(data, row['certificate'])
        entry = {'case': row['case'], 'slater': r}
        aliases = {'h6_1p6': 'fresh_h6_1p6', 'h6_1p73': 'fresh_h6_1p73'}
        name = aliases.get(row['case'], row['case'])
        if name in ('h6', 'h8', 'fresh_h6_1p6', 'fresh_h6_1p73'):
            lc = read(g.OUT/'lifted_upper'/f'{name}.json'); entry['response_lift'] = lift.check(data, lc)
            entry['lift_gain_over_HF_mHa'] = float((F(lc['hf_energy_Ha'])-F(lc['upper_Ha']))*1000)
        rows.append(entry)
    return rows


def ch2_control(uppers):
    data = read(ch2.OUT/'fixture.json'); states = []
    for spin in (0, 1):
        receipt = ch2.check(data, read(ch2.OUT/f'spin_{spin}_certificate.json'))
        U = F(next(r for r in uppers if r['case'] == f'ch2_s{spin}')['slater']['upper_Ha'])
        L = F(receipt['lower_Ha'])
        if U < L: raise ValueError('CH2 compact upper contradicts model lower')
        states.append({'spin': spin, 'compact_upper_Ha': str(U), 'inherited_lower_Ha': str(L),
                       'new_width_mHa': float((U-L)*1000), 'inherited_spin_certificate': receipt})
    Ls, Lt = [F(s['inherited_lower_Ha']) for s in states]
    Us, Ut = [F(s['compact_upper_Ha']) for s in states]
    return {'states': states, 'gap_convention': 'E_triplet-E_singlet',
            'gap_lower_Ha': str(Lt-Us), 'gap_upper_Ha': str(Ut-Ls),
            'ordering_resolved': Ut-Ls < 0 or Lt-Us > 0,
            'scope': 'Compact common-orbital HF upper controls with inherited enumerated model lower proofs; no experimental-accuracy claim.'}


def run(build=False):
    start = time.monotonic()
    closed = []
    for name in ('h6', 'h8'):
        row = reference_closed(name, build); closed.append(row)
        print(json.dumps({'stage': name, 'actions': row['terminal_aware_response']['H_actions_per_K_vector'],
                          'width_mHa': row['width_mHa']}), flush=True)
    transferred = transferred_closure(build); fresh = fresh_full_interval(); uppers = compact_uppers()
    physical = ch2_control(uppers)
    forbidden = [x for x in ('numpy', 'scipy', 'cvxpy', 'pyscf') if x in sys.modules]
    if forbidden: raise AssertionError('Numerical package loaded on exact integrated replay')
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    result = {'reference_closed': closed, 'transferred_closure': transferred, 'fresh_full_interval': fresh,
              'compact_uppers': uppers, 'ch2_control': physical, 'wall_seconds': time.monotonic()-start,
              'peak_RSS_bytes': peak, 'numerical_packages_loaded': forbidden,
              'campaign_complete_cheap_terminal_goal_achieved': False}
    if build: (g.OUT/'integrated.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'exact_replay_passed': True, 'seconds': result['wall_seconds'], 'peak_RSS_bytes': peak,
                      'fresh_width_mHa': fresh['width_mHa'], 'ch2_ordering_resolved': physical['ordering_resolved']}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--build', action='store_true')
    run(parser.parse_args().build)
