"""Independent exact family-limit replay matched to a range-two energy proof."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_range_two_family_limit import replay


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    folder = parser.parse_args().directory.resolve()
    started = time.monotonic()
    proposal_path = folder/'diagonal_family_limit_proposal.json'
    source = json.loads(proposal_path.read_text())
    if source['kind'] not in ('joint_diagonal_range2_family_limit_proposal_v1','joint_diagonal_range2_free_profile_proposal_v1','joint_diagonal_quadratic_charge_family_proposal_v1','joint_diagonal_charge_square_pair_family_proposal_v1','joint_diagonal_full_charge_indicator_family_proposal_v1','joint_diagonal_full_signed_charge_family_proposal_v1','joint_hopping_telescope_family_proposal_v1','joint_spin_telescope_family_proposal_v1','joint_spectator_hopping_family_proposal_v1','joint_pair_transfer_family_proposal_v1','joint_two_spectator_family_proposal_v1','joint_three_spectator_family_proposal_v1','joint_coherent_projector_family_proposal_v1','joint_spin_word_family_proposal_v1','joint_pure_coherence_family_proposal_v1'):
        raise ValueError('Explicit range-two proposal required')
    pure_coherence=source['kind']=='joint_pure_coherence_family_proposal_v1'
    full_spin_word=pure_coherence or source['kind']=='joint_spin_word_family_proposal_v1'
    coherent_projector=full_spin_word or source['kind']=='joint_coherent_projector_family_proposal_v1'
    three_spectator=coherent_projector or source['kind']=='joint_three_spectator_family_proposal_v1'
    two_spectator=three_spectator or source['kind']=='joint_two_spectator_family_proposal_v1'
    pair_transfer=two_spectator or source['kind']=='joint_pair_transfer_family_proposal_v1'
    spectator=pair_transfer or source['kind']=='joint_spectator_hopping_family_proposal_v1'
    spin=spectator or source['kind']=='joint_spin_telescope_family_proposal_v1'
    coherent=spin or source['kind']=='joint_hopping_telescope_family_proposal_v1'
    full_signed=coherent or source['kind']=='joint_diagonal_full_signed_charge_family_proposal_v1'
    full_indicators=full_signed or source['kind']=='joint_diagonal_full_charge_indicator_family_proposal_v1'
    square_pairs=full_indicators or source['kind']=='joint_diagonal_charge_square_pair_family_proposal_v1'
    full_quadratic=square_pairs or source['kind']=='joint_diagonal_quadratic_charge_family_proposal_v1'
    free_profile=full_quadratic or source['kind']=='joint_diagonal_range2_free_profile_proposal_v1'
    source['kind'] = 'joint_diagonal_range2_family_limit_v2' if free_profile else 'joint_diagonal_range2_family_limit_v1'
    if full_quadratic:source['kind']='joint_diagonal_range2_family_limit_v3'
    if square_pairs:source['kind']='joint_diagonal_range2_family_limit_v4'
    if full_indicators:source['kind']='joint_diagonal_range2_family_limit_v5'
    if full_signed:source['kind']='joint_diagonal_range2_family_limit_v6'
    if coherent:source['kind']='joint_diagonal_range2_family_limit_v7'
    if spin:source['kind']='joint_diagonal_range2_family_limit_v8'
    if spectator:source['kind']='joint_diagonal_range2_family_limit_v9'
    if pair_transfer:source['kind']='joint_diagonal_range2_family_limit_v10'
    if two_spectator:source['kind']='joint_diagonal_range2_family_limit_v11'
    if three_spectator:source['kind']='joint_diagonal_range2_family_limit_v12'
    if coherent_projector:source['kind']='joint_diagonal_range2_family_limit_v13'
    if full_spin_word:source['kind']='joint_diagonal_range2_family_limit_v14'
    if pure_coherence:source['kind']='joint_diagonal_range2_family_limit_v15'
    source['scope'] = 'Exact candidate input for free range-two profiles. Acceptance requires independent physical replay.' if free_profile else 'Exact candidate input for fixed range-two profile. Acceptance requires independent physical replay.'
    certificate_path = folder/'range_two_family_limit_certificate.json'
    certificate_path.write_text(json.dumps(source, indent=2)+'\n')
    result = replay(source)
    energy_path = folder/'profile_joint_r1_2_certificate.json'
    energy_receipt_path = folder/'range_two_replay.json'
    energy = json.loads(energy_path.read_text())
    if ('pure_coherence' in energy or energy.get('kind') == 'hubbard_projector_extension_v19') and not pure_coherence:
        raise ValueError('Pure-coherence energy requires a matching enlarged family cap')
    if ('spin_word_telescope' in energy or energy.get('kind') == 'hubbard_projector_extension_v18') and not full_spin_word:
        raise ValueError('Full spin-word energy requires a matching enlarged family cap')
    if ('coherent_projector' in energy or energy.get('kind') == 'hubbard_projector_extension_v17') and not coherent_projector:
        raise ValueError('Coherent-projector energy requires a matching enlarged family cap')
    if ('three_spectator_hopping' in energy or energy.get('kind') == 'hubbard_projector_extension_v16') and not three_spectator:
        raise ValueError('Three-spectator energy requires a matching enlarged family cap')
    if ('two_spectator_hopping' in energy or energy.get('kind') == 'hubbard_projector_extension_v15') and not two_spectator:
        raise ValueError('Two-spectator energy requires a matching enlarged family cap')
    if ('pair_transfer' in energy or energy.get('kind') == 'hubbard_projector_extension_v14') and not pair_transfer:
        raise ValueError('Pair-transfer energy requires a matching enlarged family cap')
    if energy.get('quadratic_charge_telescope') and not full_quadratic:
        raise ValueError('Energy quadratic telescope requires the full quadratic family cap')
    if energy.get('charge_square_pair_telescope') and not square_pairs:
        raise ValueError('Energy charge-square pairs require the matching family cap')
    if energy.get('higher_charge_indicator_telescope') and not full_indicators:
        raise ValueError('Energy higher charge indicators require the full indicator family cap')
    if energy.get('signed_charge_telescope') and not full_signed:
        raise ValueError('Energy signed-charge patterns require the full signed family cap')
    if energy.get('hopping_telescope') and not coherent:
        raise ValueError('Energy hopping telescope requires the matching family cap')
    if energy.get('spin_telescope') and not spin:
        raise ValueError('Energy spin telescope requires the matching family cap')
    if energy.get('spectator_hopping') and not spectator:
        raise ValueError('Energy spectator hopping requires the matching family cap')
    accepted = json.loads(energy_receipt_path.read_text())
    if not accepted['accepted'] or accepted['source_sha256'][str(energy_path.relative_to(ROOT))] != hashlib.sha256(energy_path.read_bytes()).hexdigest():
        raise ValueError('Energy proof does not match the current certificate')
    if set(energy['target']) != {'U', 't', 'V', 'W'} or any(F(energy['target'][k]) != v for k, v in [('U', F(4)), ('t', F(1)), ('V', F(1, 2)), ('W', F(source['W']))]):
        raise ValueError('Fixed physical target differs')
    if energy['vector'] != source['half_vector'] or energy['joint']['vector'] != source['charged_vector'] or F(energy['joint']['ratio']) != F(source['ratio']):
        raise ValueError('Projector sources or ratio differ')
    if F(energy['projector_sum_ceiling'])/energy['windows'] != F(source['theta_half']) or F(energy['joint']['projector_sum_ceiling'])/energy['joint']['windows'] != F(source['theta_joint']):
        raise ValueError('Fixed projector ceilings differ')
    if list(map(F, energy['local_window']['range_two_density_profile'])) != list(map(F, result['range_two_density_profile'])):
        raise ValueError('Range-two profile differs from the fixed dual profile')
    diagonal = {int(s): F(v) for s, v in energy['telescoping_diagonal'].items()}
    rebuilt = {}
    for shape in source['diagonal_shapes']:
        shape = {int(s): F(v) for s, v in shape.items()}
        if set(shape) & set(rebuilt):
            raise ValueError('Comparison requires disjoint supplied shapes')
        ratios = {diagonal.get(s, F(0))/v for s, v in shape.items() if v}
        if len(ratios) != 1:
            raise ValueError('Energy diagonal is outside the fixed span')
        coefficient = ratios.pop()
        rebuilt.update({s: coefficient*v for s, v in shape.items()})
    if {s: v for s, v in rebuilt.items() if v} != diagonal:
        raise ValueError('Energy has uncovered diagonal entries')
    lower = F(accepted['lower_replay']['periodic_lower_density'])
    upper = F(result['periodic_family_upper'])
    if upper < lower:
        raise ValueError('Fixed-family cap contradicts the accepted energy lower')
    files = {Path(__file__).resolve(), proposal_path, certificate_path, energy_path, energy_receipt_path}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(path).resolve())
    receipt = {'accepted': True, 'free_range_two_profile':free_profile, 'full_quadratic_charge':full_quadratic, 'charge_square_pairs':square_pairs, 'full_charge_indicators':full_indicators, 'full_signed_charge':full_signed, 'hopping_telescope':coherent, 'spin_telescope':spin, 'spectator_hopping':spectator, 'pair_transfer':pair_transfer, 'two_spectator_hopping':two_spectator, 'three_spectator_hopping':three_spectator, 'coherent_projector':coherent_projector, 'full_spin_word':full_spin_word, 'pure_coherence':pure_coherence, 'family_replay': result, 'accepted_periodic_lower': str(lower),
               'periodic_family_upper': str(upper), 'family_gap': str(upper-lower),
               'family_gap_float': float(upper-lower), 'seconds': time.monotonic()-started,
               'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
               'scope': 'Exact fixed-profile range-two family cap, matched to accepted energy sources/ceilings/ratio/shape span. Not a physical energy upper or a cap over other range-two profiles.'}
    if free_profile:receipt['scope']='Exact free range-two profile family cap, matched to accepted energy sources/ceilings/ratio/shape span. All reflected mean-correct range-two profiles covered by the additional zero moment. Not a physical energy upper or a cap over different supports.'
    if full_quadratic:receipt['scope'] += ' All six reflection-odd quadratic charge telescopes are covered by exact zero moments.'
    if square_pairs:receipt['scope'] += ' All four reflection-odd pairs of charge squares are covered by exact zero moments.'
    if full_indicators:receipt['scope'] += ' All twelve reflection-odd functions of five binary empty/double indicators are covered by exact zero moments.'
    if full_signed:receipt['scope'] += ' All 52 PH-even reflection-odd functions of five signed charges are covered.'
    if coherent:receipt['scope'] += ' The range-three hopping telescope is covered by its exact zero moment.'
    if spin:receipt['scope'] += ' All four spin-dot telescopes are covered by exact zero moments.'
    if spectator:receipt['scope'] += ' All fourteen spectator hopping telescopes are covered by exact zero moments.'
    if pair_transfer:receipt['scope'] += ' All four pair-transfer telescopes are covered by exact zero moments.'
    if two_spectator:receipt['scope'] += ' All thirty two-spectator charge-hopping telescopes are covered by exact zero moments.'
    if three_spectator:receipt['scope'] += ' All eighteen three-spectator charge-hopping telescopes are covered by exact zero moments.'
    if coherent_projector:receipt['scope']+=' Both coherent-projector overlap telescopes are covered by exact zero moments.'
    if pure_coherence:receipt['scope']+=' Both pure-coherence moments are exactly zero.'
    if full_spin_word:receipt['scope']+=' All120 diagonal spin-word moments are exactly zero.'
    (folder/'range_two_family_limit_replay.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({key: value for key, value in receipt.items() if key not in ('family_replay', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
