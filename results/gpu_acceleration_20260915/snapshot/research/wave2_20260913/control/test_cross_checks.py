"""Independent small controls for the newly introduced mathematical paths."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from experiments.marginal_transfer_verify import apply_word
from research.all_angles_20260913.frontier.exact_h4_oracle import ldl_psd
from research.wave2_20260913.tensor_network.number_conserving_mps import amplitudes, unpack, term_expect

ROOT = Path(__file__).resolve().parents[3]


def test_psd_rejects_positive_diagonal_indefinite_matrix():
    assert not ldl_psd([[F(1), F(2)], [F(2), F(1)]])[0]
    assert not ldl_psd([[F(0), F(1)], [F(1), F(0)]])[0]
    assert ldl_psd([[F(0), F(0)], [F(0), F(1)]])[0]


def test_car_transfer_matches_independent_word_action():
    data = json.loads((ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json').read_text())
    modes, particles, degeneracy = data['modes'], data['particles'], 2
    _, count = unpack(np.zeros(10000), modes, particles, degeneracy)
    parameters = np.random.default_rng(721).normal(size=count)
    vector = amplitudes(parameters, modes, particles, degeneracy)
    norm = float(vector @ vector)
    assert norm > 0
    assert np.isclose(term_expect(parameters, (), modes, particles, degeneracy), norm, rtol=1e-12)
    assert all(vector[s] == 0 for s in range(len(vector)) if s.bit_count() != particles)
    for row in data['hamiltonian']:
        word = tuple(tuple(letter) for letter in row['word'])
        exact_action_numeric_vector = 0.0
        for state, amplitude in enumerate(vector):
            if amplitude:
                target = apply_word(word, state)
                if target:
                    other, sign = target
                    exact_action_numeric_vector += amplitude*sign*vector[other]
        transfer = term_expect(parameters, word, modes, particles, degeneracy)
        assert abs(transfer-exact_action_numeric_vector)/norm < 1e-12


def test_reported_relative_width_is_endpoint_difference():
    data = json.loads((ROOT/'results/wave2_20260913/theory/parent_test.json').read_text())
    for row in data['results']:
        assert F(data['E_ref'])-F(row['lower']) == F(row['width']) == F(row['epsilon'])
