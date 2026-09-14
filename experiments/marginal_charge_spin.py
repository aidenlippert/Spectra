"""Joint charge-sector polynomials with fixed-spin coefficient bounds.

Metric kernels are computed once per charge pattern and charge-change vector,
not once per spin occupation. Charge patterns are still explicitly enumerated.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_charge_polynomial import add, multiply, scale


def cardinality_lower(poly, variables, occupied):
    """Exact lower inequality from fixed counts of active degree-d monomials."""
    if type(variables) is not int or variables < 0 or type(occupied) is not int or not 0 <= occupied <= variables:
        raise ValueError('Valid variable and occupation counts required')
    groups = {}
    for support, value in poly.items():
        if type(support) is not int or support < 0 or support >= 1 << variables:
            raise ValueError('Invalid spin monomial')
        if value: groups.setdefault(support.bit_count(), []).append(value)
    lower = F(0)
    for degree, values in groups.items():
        if degree > occupied: continue
        count = comb(occupied, degree); total = comb(variables, degree)
        negatives = sorted(x for x in values if x < 0)
        positives = sorted(x for x in values if x > 0)
        used = min(count, len(negatives))
        lower += sum(negatives[:used], F(0))
        rest = max(0, count-used-(total-len(values)))
        lower += sum(positives[:rest], F(0))
    return lower


class ChargeSpin:
    def __init__(self, data):
        self.coherent = CoherentCharge(data)
        self.signs = {}
        for key, amplitude in self.coherent.groups.items():
            c, a = key
            closed = self.coherent.close(c | a, a)
            if closed is None: continue
            lo, hi = self.coherent.amplitude_range(amplitude, *closed)
            self.signs[key] = 1 if lo >= 0 else (-1 if hi <= 0 else 0)

    def patterns(self, max_patterns=10000):
        m = self.coherent.sites
        count = sum(comb(m, d)*comb(m-d, d) for d in range(1, m//2+1))
        if type(max_patterns) is not int or not 1 <= max_patterns <= 100000 or count > max_patterns:
            raise ValueError('Charge-pattern budget exhausted or invalid')
        for d in range(1, m//2+1):
            for plus in combinations(range(m), d):
                for minus in combinations([i for i in range(m) if i not in plus], d):
                    yield tuple(1 if i in plus else (-1 if i in minus else 0) for i in range(m))

    def compile(self, charges):
        o = self.coherent
        if (type(charges) not in (list, tuple) or len(charges) != o.sites or
                any(type(q) is not int or q not in (-1, 0, 1) for q in charges) or
                sum(charges) or 1 not in charges):
            raise ValueError('Complete neutral ionic charge pattern required')
        singles = [i for i, q in enumerate(charges) if q == 0]
        r = len(singles); positions = {i: j for j, i in enumerate(singles)}
        one = {0: F(1)}
        n = []
        for i, q in enumerate(charges):
            if q: n.extend([one if q == 1 else {}]*2)
            else:
                x = {1 << positions[i]: F(1)}
                n.extend([x, add(one, scale(x, -1))])
        row = {}
        for support, value in o.oracle.diagonal.items():
            term = one
            for i in range(o.modes):
                if support & (1 << i): term = multiply(term, n[i])
            row = add(row, scale(term, value))
        kernels = {}; groups = metric_evaluations = 0
        for (c, a), amplitude in o.groups.items():
            event = one
            for i in range(o.modes):
                if c & (1 << i): event = multiply(event, add(one, scale(n[i], -1)))
                elif a & (1 << i): event = multiply(event, n[i])
            if not event: continue
            delta = tuple(((c >> (2*i)) & 3).bit_count()-((a >> (2*i)) & 3).bit_count()
                          for i in range(o.sites))
            changed = tuple(q+d for q, d in zip(charges, delta))
            if not any(changed): continue
            amp = {}
            for support, value in amplitude.items():
                amp = add(amp, scale(n[support.bit_length()-1] if support else one, value))
            joint = multiply(event, amp)
            if not joint: continue
            sign = self.signs.get((c, a), 0)
            if not sign:
                raise ValueError('Grouped amplitude changes sign; this compiler requires a certified sign')
            if delta not in kernels:
                ratio = F(1)
                for label, factor in o.local_factors:
                    if not any(delta[i] for i, _ in label): continue
                    before = after = 1
                    for i, power in label:
                        before *= charges[i]**power; after *= changed[i]**power
                    ratio *= factor**(after-before); metric_evaluations += 1
                kernels[delta] = ratio
            row = add(row, scale(joint, -sign*kernels[delta])); groups += 1
        # Monomials with degree>r/2 vanish on the admitted spin slice.
        row = {mask: value for mask, value in row.items() if mask.bit_count() <= r//2}
        return row, {'single_sites': r, 'spin_assignments': comb(r, r//2),
                     'transition_groups': groups, 'metric_kernels': len(kernels),
                     'metric_scalar_evaluations': metric_evaluations, 'polynomial_terms': len(row),
                     'maximum_degree': max((mask.bit_count() for mask in row), default=0),
                     'spin_endpoint_evaluations': 0}


def quadratic_slice_gate(poly, variables, target):
    """Exact projected PSD certificate on sum(z)=0, z_i^2=1.

    The optional four-site diagonal gauge is deterministic. It diagonalizes
    the balanced subspace, so the gate is exact there for even quadratic rows.
    For larger quadratic rows the zero gauge gives only a sufficient bound.
    """
    from experiments.marginal_schur_transfer import ldl_pivots
    r = variables
    if type(r) is not int or r < 0 or r % 2 or any(type(mask) is not int or not 0 <= mask < 1 << r for mask in poly):
        raise ValueError('Even balanced spin slice with valid monomials required')
    if r < 2 or any(mask.bit_count() > 2 for mask in poly): return None
    c = poly.get(0, F(0))
    linear = [F(0)]*r; matrix = [[F(0) for _ in range(r)] for _ in range(r)]
    for mask, coefficient in poly.items():
        sites = [i for i in range(r) if mask & (1 << i)]
        if len(sites) == 1:
            i, = sites; c += coefficient/2; linear[i] += coefficient/2
        elif len(sites) == 2:
            i,j = sites; c += coefficient/4
            linear[i] += coefficient/4; linear[j] += coefficient/4
            matrix[i][j] += coefficient/8; matrix[j][i] += coefficient/8
    if any(value != linear[0] for value in linear): return None
    diagonal = [F(0)]*r
    if r == 4:
        basis = [(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)]
        for a,b in combinations(basis,2):
            cross = sum((a[i]*matrix[i][j]*b[j] for i in range(r) for j in range(r)),F(0))
            for i in range(r): diagonal[i] -= cross*a[i]*b[i]/4
    assert sum(diagonal)==0
    threshold = (F(target)-c)/r
    shifted = [[matrix[i][j]+((diagonal[i]-threshold) if i==j else 0) for j in range(r)] for i in range(r)]
    # B columns are e_i-e_last; B spans the full balanced subspace.
    projected = [[shifted[i][j]-shifted[i][-1]-shifted[-1][j]+shifted[-1][-1]
                  for j in range(r-1)] for i in range(r-1)]
    pivots = ldl_pivots(projected)
    if pivots is None: return None
    return {'matrix_dimension':r-1,'pivots':[str(x) for x in pivots],
            'diagonal_gauge':[str(x) for x in diagonal],'spectral_threshold':str(threshold)}


def replay(certificate):
    if certificate.get('kind') != 'valence_charge_spin_v1':
        raise ValueError('Unsupported charge-spin certificate')
    if type(certificate.get('target_lower')) is not str:
        raise ValueError('Exact requested lower bound required')
    target = F(certificate['target_lower']); oracle=ChargeSpin(certificate)
    max_patterns=certificate.get('max_charge_patterns')
    stats={'charge_patterns':0,'covered_Q_configurations':0,'cardinality_certificates':0,
           'quadratic_psd_certificates':0,'maximum_psd_dimension':0,'metric_scalar_evaluations':0,
           'metric_kernels':0,'transition_groups':0,'polynomial_terms':0,'maximum_polynomial_degree':0,
           'spin_endpoint_evaluations':0,'determinant_actions':0}
    for charges in oracle.patterns(max_patterns):
        poly,cost=oracle.compile(charges)
        r=cost['single_sites']; lower=cardinality_lower(poly,r,r//2)
        if lower>=target: stats['cardinality_certificates']+=1
        else:
            gate=quadratic_slice_gate(poly,r,target)
            if gate is None: raise ValueError('Charge-sector spin polynomial misses requested bound')
            stats['quadratic_psd_certificates']+=1
            stats['maximum_psd_dimension']=max(stats['maximum_psd_dimension'],gate['matrix_dimension'])
        stats['charge_patterns']+=1; stats['covered_Q_configurations']+=cost['spin_assignments']
        for key in ['metric_scalar_evaluations','metric_kernels','transition_groups','polynomial_terms']:
            stats[key]+=cost[key]
        stats['maximum_polynomial_degree']=max(stats['maximum_polynomial_degree'],cost['maximum_degree'])
    o=oracle.coherent; expected=o.oracle.sector_dimension-comb(o.sites,o.target)
    if stats['covered_Q_configurations']!=expected: raise ValueError('Incomplete ionic charge coverage')
    stats.update(complement_lower=str(target),source_hamiltonian_terms=len(o.oracle.h),
                 coherent_transition_groups=len(o.groups),fixed_sign_groups=sum(bool(x) for x in oracle.signs.values()),
                 scope='Complete explicit ionic charge-pattern coverage with exact grouped row polynomials on the balanced spin slice. Metric kernels are shared across spin occupations; coefficient-count or projected PSD inequalities certify each whole slice. No spin endpoint evaluation or determinant action. Charge-pattern count remains combinatorial; the four-spin PSD gauge is special and general scaling is not established.')
    return stats


def build(source, output, target, max_patterns=10000):
    import json
    from pathlib import Path
    source,output=Path(source),Path(output)
    if output.exists(): raise ValueError('Preserve previous charge-spin proof')
    data=json.loads(source.read_text())
    certificate={k:data[k] for k in ['modes','particles','hamiltonian','metric_rule']}
    certificate.update(kind='valence_charge_spin_v1',target_lower=str(F(target)),max_charge_patterns=max_patterns)
    receipt=replay(certificate)
    output.mkdir(parents=True)
    for name,obj in [('certificate',certificate),('receipt',receipt)]:
        (output/(name+'.json')).write_text(json.dumps(obj,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source'); p.add_argument('--output'); p.add_argument('--verify')
    p.add_argument('--target',default='-6.264'); p.add_argument('--max-charge-patterns',type=int,default=10000)
    args=p.parse_args()
    print(json.dumps(replay(json.loads(Path(args.verify).read_text())) if args.verify else
                     build(args.source,args.output,args.target,args.max_charge_patterns),indent=2))
