"""Maintain a bounded CH2 model-error branch using preserved numerical controls."""
from fractions import Fraction as F
import hashlib
import json
from research.interacting_scaling_20260915.budget import ROOT, OUT, dump


def run():
    sources = {}
    def read(path):
        sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    conversion = 627.509474
    old = read(ROOT/'results/correlated_pair_20260913/ch2/interval.json')
    states = {row['spin']: row for row in old['states']}
    lower = F(states[0]['lower_Ha'])-F(states[1]['upper_Ha'])
    upper = F(states[0]['upper_Ha'])-F(states[1]['lower_Ha'])
    if lower != F(old['gap_lower_Ha']) or upper != F(old['gap_upper_Ha']):
        raise ValueError('Preserved spin-gap subtraction does not reproduce')
    rows = []
    base = ROOT/'results/transfer_solver_20260915/ch2_model_study'
    for basis in ('cc-pvdz', 'cc-pvtz'):
        controls = [read(base/f'{basis}_S{spin}'/'result.json') for spin in (0, 1)]
        cc = [read(base/f'{basis}_S{spin}'/'cc_control.json') for spin in (0, 1)]
        if any(d['total_spin'] != spin or d['basis'] != basis or not d['spin_checked_at_every_evaluation']
               for spin, d in enumerate(controls)):
            raise ValueError('State/model identification did not reproduce')
        for method, key, values in [('CASSCF(6,6)', 'CASSCF_total_Ha', controls),
                                    ('SC-NEVPT2', 'SC_NEVPT2_total_Ha', controls),
                                    ('CCSD(T)', 'CCSD_T_total_Ha', cc)]:
            rows.append({'basis': basis, 'method': method,
                'electronic_gap_kcal_mol': (values[0][key]-values[1][key])*conversion,
                'geometries': [d['geometry_Angstrom'] for d in controls],
                'geometries_optimized_for': 'CASSCF(6,6) in this basis',
                'solver_interval_certified': False, 'ZPE_computed': False})
    target = {'observable': 'CH2 E(1A1)-E(3B1); distinguish relaxed electronic Te from the vibrational ground-state splitting',
        'preserved_fixed_geometry_model': {'gap_lower_Ha': str(lower), 'gap_upper_Ha': str(upper),
            'display_gap_kcal_mol': [float(lower)*conversion, float(upper)*conversion],
            'solver_width_kcal_mol': float(upper-lower)*conversion,
            'lower_dependency': 'Enumerated spin-projection matrices; retained as a model control',
            'physical_error_interval': None},
        'numerical_model_ladder': rows,
        'experimental_context': {'inferred_Te_kcal_mol': [8.7, .5], 'observed_splitting_kcal_mol': [9., .09],
            'source': 'https://doi.org/10.1063/1.449746', 'checked_primary_abstract': '2026-09-15'},
        'missing_physical_contributions': ['Controlled complete-basis and active-space limits',
            'Geometry errors for the post-CASSCF methods', 'Zero-point and anharmonic nuclear corrections',
            'Remaining electronic/relativistic contributions at the chosen comparison tolerance'],
        'no_combined_certified_physical_error_bar': True,
        'basis_comparison_also_changes_optimized_geometry': True,
        'no_new_chemistry_or_prospective_predictions_in_this_small_branch': True,
        'source_sha256': sources}
    dump(OUT/'physical_model_ladder.json', target)


if __name__ == '__main__': run()
