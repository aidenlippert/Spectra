"""Final outcome supplement; run before the frozen report seals this pass."""
from fractions import Fraction as F
import json
from pathlib import Path
from research.transfer_solver_20260915.budget import ROOT, OUT, dump


def read(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def run():
    paths = [OUT/'cases'/n for n in ('h4_control', 'h8_cold', 'h6_asymmetric', 'water_asymmetric')]
    attempts = [read(p) for p in sorted((OUT/'runs').glob('*.json'))]
    byname = {x['name']: x for x in attempts}
    rows, results = [], []
    for p in paths:
        interval = read(p/'interval.json')
        if p.name == 'h10_rotated_upper':
            stages = [x for x in attempts if x['name'].startswith(('h10_rotated_lower_', 'h10_rotated_'))]
            # The prefixes overlap intentionally; each receipt is selected once.
            stages.append(byname['h10_size_generate'])
        else:
            stages = [byname[p.name+'_'+n] for n in ('generate', 'state', 'upper', 'nonsinglet', 'prepare', 'stage1', 'stage2', 'replay')]
        cost = sum(x['wall_seconds'] for x in stages)
        peak = max(x['peak_child_RSS_bytes'] for x in stages)
        width = f"{interval['width_mHa']:.9f}" if interval else 'No complete interval'
        forecast = read(p/'prepared/forecast.json', {})
        rows.append(f"| {p.name} | {width} | {cost:.3f} | {peak/1e6:.1f} | {forecast.get('retained_Gram_entries', '—')} |")
        results.append({'case': p.name, 'interval': interval, 'stage_process_seconds': cost,
            'peak_RSS_bytes': peak, 'forecast': forecast, 'execution': read(p/'execution.json')})
    upper = read(OUT/'rotated_h10/original_upper.json')
    nvidia = read(ROOT/'results/nvidia_followup_20260915/summary.json', {})
    reference = read(OUT/'cases/h10_size/fci_reference.json')
    h10 = read(OUT/'adaptive/h10_rotated_upper/interval.json')
    h10_lower = read(OUT/'adaptive/h10_rotated_upper/lower.json')
    upper_error = 1000*(float(F(upper['upper_Ha']))-reference['numerical_electronic_energy_Ha'])
    h10_text = 'The full H10 interval was not completed; the accepted upper remains a separate achieved component.'
    if h10:
        lower_error = 1000*(reference['numerical_electronic_energy_Ha']-float(F(h10['lower_Ha'])))
        h10_text = f"The exact full H10 width is **{h10['width_mHa']:.9f} mHa**. The numerical reference lies {lower_error:.9f} mHa above the lower and {upper_error:.9f} mHa below the upper."
        if h10_lower:
            h10_text += f" The singlet residual allowance is {1000*float(F(h10_lower['singlet']['residual_l1'])):.9f} mHa."
    h10_text += ' This is the HF-guided lower control, not the later correlated-moment result.'
    correlated = nvidia.get('best_H10')
    correlated_text = 'A later correlated-moment GPU attempt is recorded separately; no new exact outcome is asserted here without its complete replay receipt.'
    if correlated:
        correlated_text = f"The correlated-moment follow-up obtained a complete exact width of **{correlated['width_mHa']:.9f} mHa** on the original H10 Hamiltonian. The exact upper is the same one described above. Its source, optimizer, library backends, independent replay, complete dependency accounting, and unsuccessful attempts are in the [NVIDIA follow-up report]({ROOT}/research/nvidia_followup_20260915/REPORT.md)."
    recovery = read(OUT/'controls/h8_restricted64/recovered/interval.json')
    recovery_text = ('The saved restricted-family point was not independently accepted.' if not recovery else
        f"The saved restricted64 iterate was independently accepted at **{recovery['width_mHa']:.9f} mHa**. This is a valid attained bound, not an optimum or a family ceiling.")
    global_result = read(OUT/'fragments_global_charge/summary.json')
    global_text = 'The attempted global-charge extension did not complete; inspect its execution receipt.'
    if global_result:
        global_text = rf"""The new charge-sector control removes that local-charge restriction for these
H4 fragments. It independently verifies the other eight local charge sectors,
enumerating **{global_result['extra_local_determinant_labels_enumerated']} local determinant labels**, and establishes
\(H_i\ge L_4+\mu(N_i-4)\), with \(\mu\approx {float(F(global_result['chemical_potential_Ha'])):.9f}\) Ha.
On the complete total-number sector \(\sum_i N_i=4k\), summing gives
\(\sum_iH_i\ge kL_4\), even when local charge redistributes.

For 2, 4 and 8 fragments the exact interval widths remain 0.003089495,
0.006178991 and 0.012357982 mHa. Neither composite determinants nor charge
assignments were enumerated. Missing charge sectors, an invalid supporting
slope and interfragment coupling are rejected. This establishes this composite
certificate rule; it does not test a monolithic optimizer's size consistency."""
    physical_rows, physical = [], []
    for b in ('cc-pvdz', 'cc-pvtz'):
        s, t = [read(OUT/f'ch2_model_study/{b}_S{k}/result.json') for k in (0, 1)]
        cs, ct = [read(OUT/f'ch2_model_study/{b}_S{k}/cc_control.json') for k in (0, 1)]
        if s and t:
            casscf = (s['CASSCF_total_Ha']-t['CASSCF_total_Ha'])*627.509474
            nevpt = (s['SC_NEVPT2_total_Ha']-t['SC_NEVPT2_total_Ha'])*627.509474
            cc = (cs['CCSD_T_total_Ha']-ct['CCSD_T_total_Ha'])*627.509474 if cs and ct else None
            cc_text = f'{cc:.5f}' if cc is not None else 'not completed'
            physical_rows.append(f'| {b} | {casscf:.5f} | {nevpt:.5f} | {cc_text} |')
            physical.append({'basis': b, 'CASSCF_gap_kcal_mol': casscf, 'SC_NEVPT2_gap_kcal_mol': nevpt, 'CCSD_T_gap_kcal_mol': cc})
    summary = {'cases': results, 'H10_rotated_upper_reference_excess_mHa': upper_error,
        'H10_upper_rotation_allowance_mHa': 1000*float(F(upper['rotation_rounding_allowance_Ha'])),
        'restricted64_recovered_interval': recovery, 'global_fragment_control': global_result,
        'physical_method_comparison': physical,
        'NVIDIA_followup': nvidia,
        'all_instrumented_attempts': len(attempts),
        'all_instrumented_process_wall_seconds': sum(x['wall_seconds'] for x in attempts),
        'all_open_research_questions_solved': False}
    dump(OUT/'final_outcomes.json', summary)
    write_document(summary, rows, h10_text, correlated_text, recovery_text, global_text, physical_rows)


def write_document(summary, rows, h10_text, correlated_text, recovery_text, global_text, physical_rows):
    text = rf'''# Verified outcomes and remaining research questions

The preserved H8 milestone remains valid. This pass adds **cold construction,
successful limited transfer, explicit comparisons, and a complete H10 certificate
through the correlated-moment follow-up**. The achieved results and unresolved questions below have different
evidence requirements; a finite campaign cannot certify broad usefulness.

## Completed model certificates

| Case | Full fixed-N width (mHa) | Successful-path process seconds from integrals | Peak process MB | Singlet Gram entries |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

Times include new integrals, state construction, nonsinglet screening,
preparation, optimization and exact acceptance on the listed successful path.
The later H10 result and its orbital-rotation upper are detailed below. Failed
preceding attempts and comparison work are additionally charged in the complete ledger.
No full fixed-N determinant space is enumerated by these molecular constructors
or accepting checks. They do evaluate a single HF determinant and construct
finite-degree operator coefficient and Gram arrays.

H8 is a fresh run from integrals without old states, coefficient maps, winning
checkpoints or reference vectors. Its canonical orbital signs differ from the
original H8 only by two exactly verified paired sign flips; the coefficient
residual of that bridge is zero. The uniform dictionary adds mixed linear words
and has 343,424 Gram entries versus 335,168 in the older winning construction.
Thus this is a cold success of the declared input-derived rule, not a claim that
the old optimization has been replicated with precisely identical parameters.

The asymmetric H6 case loses the nontrivial spatial parity exploited by H8.
Water changes both the molecule and the electron/orbital counts. Both were
frozen before their outcomes were known. This establishes transfer on these
cases; it does not establish a statistical success rate over unseen chemistry.

## H10: the state bottleneck and the full lower

The original random-state H10 run timed out after 900 seconds. A separately
declared smaller random-state retry timed out after 600 seconds. Their costs
and completed-sweep checkpoints remain preserved. Neither is counted as a
success of the original frozen rule.

An adaptive local-orbital construction then produced an upper for the **same
rational H10 Hamiltonian**. It uses a bond-32 MPS with 3,904 nonzero integer tensor
entries, a freshly made alternating-spin product initializer, and no FCI input.
The generic initializer label in the imported discovery module is superseded
by its recorded `initial_MPS_path` and the actual new `neel_seed.json`.

The checked upper is about **{summary['H10_rotated_upper_reference_excess_mHa']:.9f} mHa** above the tight numerical FCI reference.
The orbital-coefficient rounding allowance is only
**{summary['H10_upper_rotation_allowance_mHa']:.9f} mHa**. Rotation construction, state
discovery and upper replay took approximately 296 seconds, plus the charged
ordinary integral generation. The numerical DMRG convergence flag was false;
the independently evaluated exact Rayleigh quotient supplies the guarantee.
The numerical penalized sweep energy is not the accepted upper.

{h10_text}

The HF-guided lower retains the original canonical-basis dictionary and exact checker.
A newly made HF product supplies proposal moments only; the upper comes from
the separately checked correlated state. This separates upper quality from
lower-family performance. The H10 Gram entry count is about 3.75 times H8's,
and its largest block is 505 versus 260. That is measured growth, not a scaling
theorem or a demonstration of large-molecule practicality.

The next attempt transported reduced density matrices through order three from
the new local-orbital MPS into the original canonical basis. That took 15.482
process seconds, enumerated no full-N determinants, passed all fixed-charge
trace checks, and reproduced the physical energy within 3.61e-9 Ha of the checked
upper. It retained the same coefficient maps and started new Gram coordinates
at zero. Map preparation, the MPS and nonsinglet proof remain charged inputs.
The planned local two-stage correlated solve was superseded by the explicitly
authorized isolated GPU experiment; it was not executed.

{correlated_text}

The complete rotation argument is in
[ORBITAL_ROTATION.md]({ROOT}/research/transfer_followup_20260915/ORBITAL_ROTATION.md).
It applies an exactly orthogonal rational transformation to one- and two-orbital
coefficient arrays and pays an exact coefficient-norm bound. It does not form
the many-body unitary or fit a compressed state from an FCI vector.

## What the representation controls settle

The unreduced magnetic H8 control attained **94.301434775 mHa**. The accepted
0.605110982 mHa certificate is an exact member of that control's feasible
family, including its ideal-multiplier spans. Thus the poor magnetic result
cannot establish that its family lacks a sufficiently good certificate.
It demonstrates a finite-budget optimization/conditioning limitation for this
representation. The membership check uses rational span elimination and binds
the already accepted source certificate by hashes.

The restricted64 run stopped at the unchanged numerical coefficient
reconstruction guard. Its last saved point was separately exported and sent
through the original complete exact checker.

{recovery_text}

An exact dual obstruction for that restricted singlet family is still absent.
The older fixed-number trace repair cannot simply be reused: the present
singlet dual requires \(y(S^2)=0\), whereas the normalized full fixed-N trace
on m spin orbitals has
\(\tau_N(S^2)=3N(m-N)/[4(m-1)]\). For H8 it equals 3.2. Mixing in that trace
would violate the singlet ideal. A valid dual repair must respect the singlet
identities and their forced null spaces. This explains a mathematical obstacle
to reusing that repair; it is not an impossibility result for the restricted cone.

## Independent validation, composition and practical cost

A separately implemented H4 bit-state Hamiltonian and rational LDL check
accept the full lower; independent integer state expansion reproduces the upper.
A deliberately raised lower is rejected. This validation explicitly enumerates
70 small-system states and is charged separately.

The first fragment check certified disjoint sums with fixed local charges.
{global_text}

The exact enumerated baselines meet the same 1.6 mHa target much faster on H4,
asymmetric H6 and water. Their actual widths, state counts and complete costs
are in the [frozen campaign report]({ROOT}/research/transfer_solver_20260915/REPORT.md).
Numerical FCI is also much faster on these small cases, but does not provide the
same two-sided exact guarantee. This pass demonstrates no general competitive
advantage. The faster spin projector was checked coefficient-for-coefficient;
it is a verified component improvement, not a full-pipeline speedup measurement.

## Physical-model comparison

All entries below use \(E_S-E_T\) in kcal/mol. CASSCF optimizes the two states'
geometries separately. SC-NEVPT2 and CCSD(T) are single points at those same
CASSCF geometries. They are conventional numerical comparisons, not rigorous
Spectra physical predictions.

| Basis | CASSCF(6,6) | SC-NEVPT2 | CCSD(T) at CASSCF geometries |
|---|---:|---:|---:|
{chr(10).join(physical_rows)}

The measured splitting is 9.00 +/- 0.09 kcal/mol; the study infers an electronic
separation of 8.7 +/- 0.5 kcal/mol after assessing zero-point effects.
[Leopold et al.](https://experts.umn.edu/en/publications/methylene-a-study-of-the-xsup3supbsub1sub-and-%C3%A3sup1supasub1sub-st/).
The new CASSCF and SC-NEVPT2 controls still exceed that electronic reference.
The old precisely solved minimal-basis gap near 25.49 kcal/mol was therefore not
a reliable prediction of this observable. Basis and method changes are
sensitivity measurements, not certified physical-error bounds.

CCSD(T) correlates all eight electrons here and is nonvariational. Its triplet
uses ROHF orbitals in unrestricted CC storage; exact spin purity of the CC state
has not been certified. No CCSD(T) geometry optimization or zero-point correction
was performed. See the [PySCF method documentation](https://pyscf.org/user/cc.html).
No prospective experimental choice or Spectra-specific predictive advantage has
been demonstrated by these retrospective controls.

## Accounting and scope

The local computation in this pass records **{summary['all_instrumented_attempts']} instrumented attempts** totaling
**{summary['all_instrumented_process_wall_seconds']:.3f} process-wall seconds**. The complete
ledger includes failures, failed searches, reference comparisons and verification.
It separates that campaign total from the cost of a successful fresh solution.
Editing, reading, web research, reporting, and preservation I/O are outside this
compute ledger; earlier R&D remains charged in its preserved prior ledgers.
The separate NVIDIA report accounts for remote preparation, installation,
compute, verification, transfer, all attempts and the owned instance lifetime.
Mixed-host component sums are distinguished from a single freshly timed complete
pipeline. The earlier GPU acceleration task's files and results are preserved.

Four initial FCI comparison exports failed on a NumPy boolean. Their partial
files were preserved. Four repaired runs wrote valid records, then failed while
printing the same boolean. A separate validation receipt confirms the complete
records and their numerical agreement; process failures were not erased.
Four initial CH2 dispatches were rejected before chemistry started because their
job names used uppercase letters. Corrected dispatches all completed. The first
chemistry job spent substantial wall time in runtime loading; that delay remains
in its measured cost. Frozen scientific sources and accepting gates were retained.

The predeclared canonical-HF state retry and generic `refinements run` helper
were not executed. They are not evidence. The actual adaptive H10 path is the
orbital-rotation construction and canonical lower described above.

Cold discovery, these transfer cases, small independent validation and the
explicit controls now have concrete receipts. Broad transfer reliability,
competitive cost on larger correlated systems, a restricted-family optimum,
monolithic optimizer size consistency, certified physical-model accuracy and
prospective experimental usefulness remain research questions. None is labeled
solved merely because a small benchmark passes.

- [Machine-readable final outcomes]({OUT}/final_outcomes.json)
- [All-attempt compute ledger]({OUT}/cost_ledger.json)
- [Exact proof derivation]({ROOT}/research/transfer_solver_20260915/DERIVATION.md)
- [Preservation and seal]({OUT}/seal.json)
'''
    path = ROOT/'research/transfer_followup_20260915/RESULTS.md'
    with path.open('x') as stream:
        stream.write(text)
    print(path)


if __name__ == '__main__':
    run()
