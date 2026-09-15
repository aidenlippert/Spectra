"""Certify the energetic cost of retaining the uncoupled H8 state at full coupling."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, dump
from research.interacting_scaling_20260915.coupling import initialize as coupled_case
from research.interacting_scaling_20260915.pipeline import execute, PREFIX

NAME = 'h8_frozen_uncoupled_state'


def initialize():
    zero = OUT/'cases/h8_adaptive_coupling_l0_base'
    source = OUT/'cases/h8_matched_direct_base'
    case = coupled_case(source, NAME, F(1), widths=[2], previous=zero)
    (case/'mps').mkdir()
    shutil.copyfile(case/'continuation_seed.json', case/'mps/state.json')
    old = json.loads((zero/'mps/state.json').read_text())
    new = json.loads((case/'mps/state.json').read_text())
    if any(old[k] != new[k] for k in ('tensors', 'bond_charges', 'denominator', 'spin_counts')):
        raise ValueError('The reference state changed')
    dump(case/'control.json', {'old_state': str(zero/'mps/state.json'),
        'old_state_sha256': hashlib.sha256((zero/'mps/state.json').read_bytes()).hexdigest(),
        'wavefunction_unchanged': True, 'only_Hamiltonian_fixture_binding_changed': True,
        'source_state_discovery_additional': True, 'not_a_fresh_state_search': True})


def compare():
    case = OUT/'cases'/NAME
    path = json.loads((OUT/'coupling/h8_adaptive_coupling/results.json').read_text())
    best = next((row.get('best') for row in path['outcomes'] if row['lambda'] == '1'), None)
    reference = Path(best['case']) if best else OUT/'cases/h8_matched_direct_level0'
    if json.loads((case/'fixture.json').read_text()) != json.loads((reference/'fixture.json').read_text()):
        raise ValueError('The two energy calculations have different Hamiltonians')
    bound = json.loads((reference/'original_interval.json').read_text())
    energy = F(json.loads((case/'upper.json').read_text())['upper_Ha'])
    eta = F(bound['rotation_allowance_each_endpoint_Ha'])
    lower = energy-eta-F(bound['upper_Ha'])
    upper = energy+eta-F(bound['lower_Ha'])
    if lower > upper: raise ValueError('Reversed energy-excess interval')
    dump(OUT/'coupling/h8_adaptive_coupling/frozen_state_control.json', {
        'comparison_ground_certificate': str(reference), 'frozen_state_case': str(case),
        'frozen_state_energy_local_Ha': str(energy),
        'original_H_frozen_state_excess_interval_Ha': [str(lower), str(upper)],
        'original_H_frozen_state_excess_interval_mHa': [float(1000*lower), float(1000*upper)],
        'rotation_allowance_on_reference_state_and_ground_interval_included': True,
        'meaning': 'Rigorous excess energy of the rotated uncoupled reference state over the fully coupled original-model ground energy',
        'not_a_timing_advantage_or_a_new_ground_certificate': True,
        'source_state_and_ground_certificate_costs_additional': True})


def run():
    return execute(NAME, [('initialize', 30, [STD, '-B', '-S', '-m', PREFIX+'response_control', 'initialize']),
        ('upper', 180, [STD, '-B', '-S', '-m', 'research.transfer_solver_20260915.actions', str(OUT/'cases'/NAME), 'upper']),
        ('compare', 30, [STD, '-B', '-S', '-m', PREFIX+'response_control', 'compare'])])


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('initialize', 'compare', 'run'))
    a = p.parse_args(); {'initialize': initialize, 'compare': compare, 'run': run}[a.action]()
