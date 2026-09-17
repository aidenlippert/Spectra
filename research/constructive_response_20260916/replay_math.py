"""Standard-library replay of accepted finite results, including corruptions."""
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
from research.constructive_response_20260916 import dimer_family
from research.constructive_response_20260916.bond_coercivity import verify_bond_all_sectors
from research.constructive_response_20260916.composition_exact import exact_composition,nested_exact,additive_error_counterexample
from research.constructive_response_20260916.metric_certificate import check_square,sharpness_counterexample


def replay(certificate,out):
    if out.exists():
        raise FileExistsError(out)
    start=time.monotonic();out.mkdir(parents=True)
    payload=json.loads(certificate.read_text());accepted=check_square(payload)
    corruptions=[]
    for name in ('zero_rho','zero_trial','wrong_reference','wrong_model'):
        broken=copy.deepcopy(payload)
        if name=='zero_rho': broken['rho_over_t']='0'
        if name=='zero_trial': broken['retained_trial']=['0']*len(broken['retained_trial'])
        if name=='wrong_reference': broken['reference_energy_over_t']='0'
        if name=='wrong_model': broken['kind']='different_model'
        try:
            check_square(broken)
        except ValueError as err:
            corruptions.append({'case':name,'rejected':True,'reason':str(err)})
        else:
            raise AssertionError('Corruption accepted: '+name)
    for t in (F(0),F(1,2),F(1),F(2)):
        for a in (0,1,2,3,4,8):
            if not verify_bond_all_sectors(F(a),t)[0]:
                raise AssertionError('Local bond PSD verification failed')
    abstract={'response':exact_composition(),'nested':nested_exact(),
              'shift_amplification':additive_error_counterexample(),
              'sharp_metric_constant':sharpness_counterexample(F(2),F(1,10),F(9999,1800000)),
              'no_uniform_bound_at_threshold':sharpness_counterexample(F(2),F(1),F(10**6))}
    dimers=dimer_family.run()
    for family in dimers['families']:
        for row in family['rows']:
            if not 0<F(row['global_Qgap_lower_over_t'])<=F(row['global_Qgap_upper_over_t']):
                raise AssertionError('Invalid product-family gap enclosure')
    (out/'square.json').write_text(json.dumps(accepted,indent=2)+'\n')
    (out/'exact_composition.json').write_text(json.dumps(abstract,default=str,indent=2)+'\n')
    (out/'dimer_family.json').write_text(json.dumps(dimers,indent=2)+'\n')
    result={'status':'accepted_scoped_exact_results_only','certificate_sha256':hashlib.sha256(certificate.read_bytes()).hexdigest(),
            'square_model':accepted,'local_bond_checks':24,'corruption_checks':corruptions,
            'global_states_constructed_for_dimer_family':0,
            'global_states_enumerated_for_square':36,
            'new_field_level_breakthrough_established':False,
            'seconds':time.monotonic()-start}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    replay(a.certificate,a.out)
