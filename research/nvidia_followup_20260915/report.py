"""Summarize measured library effects, certified outcomes and complete cloud cost."""
from datetime import datetime
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/nvidia_followup_20260915'
CASES=ROOT/'results/transfer_solver_20260915'


def read(path,default=None):
    return json.loads(path.read_text()) if path.exists() else default


def process_metrics(path):
    if not path.exists():return {}
    text=path.read_text()
    result={}
    for line in text.splitlines():
        if 'Elapsed (wall clock)' in line:
            parts=[float(v) for v in line.rsplit(': ',1)[1].split(':')]
            result['wall_seconds']=sum(v*60**i for i,v in enumerate(reversed(parts)))
    match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)
    if match:result['peak_host_RSS_bytes']=int(match[1])*1024
    return result


def run():
    if not (OUT/'termination_confirmed.json').exists():
        raise ValueError('Confirm owned-instance termination before final accounting')
    remote=OUT/'download_final/results/nvidia_followup_20260915'
    owned=read(OUT/'instance.json')
    stopped=read(OUT/'termination_confirmed.json')
    lifetime=(datetime.fromisoformat(stopped['confirmed_UTC'])-datetime.fromisoformat(owned['launch_UTC'])).total_seconds()
    cloud={'owned_instance_id':owned['id'],'launch_to_absence_seconds':lifetime,
        'price_USD_per_hour':owned['price_cents_per_hour']/100,
        'estimated_USD_including_boot_and_idle':lifetime/3600*owned['price_cents_per_hour']/100,
        'actual_provider_invoice_not_available':True,'owned_instance_terminated':True,'preexisting_instances_untouched':True}
    intervals=[]
    h10=CASES/'adaptive/h10_correlated_guide'
    for name in ('hybrid_adaptive','sparse_cudss_adaptive'):
        exact=h10/name/'exact_local_stdlib'
        interval=read(exact/'interval.json')
        if interval:
            if not read(exact/'complete_replay.json')['all_upper_and_lower_dependencies_rechecked']:
                raise ValueError('Complete independent replay missing')
            compiled=read(h10/name/'exact_compiled_retry/interval.json') or read(h10/name/'exact_compiled/interval.json')
            if compiled and any(F(compiled[key])!=F(interval[key]) for key in ('lower_Ha','upper_Ha')):
                raise ValueError('Compiled and standard-library endpoints differ')
            intervals.append({'proposal':name,'path':str(exact),'width_mHa':interval['width_mHa'],
                'lower_Ha':interval['lower_Ha'],'upper_Ha':interval['upper_Ha'],'target_met':interval['target_met'],
                'local_complete_replay':read(exact/'complete_replay.json')})
    if not intervals:
        raise ValueError('A locally rechecked complete H10 certificate is required')
    best=min(intervals,key=lambda x:F(x['upper_Ha'])-F(x['lower_Ha']))
    comparisons=[]
    job_rows=read(remote/'accelerated_jobs.json',[])
    by_job={row['name']:row for row in job_rows}
    for label,path in [('H8',CASES/'cases/h8_cold'),('H10',h10)]:
        rows=[]
        for tag in ('dense_200','sparse_cpu_200','sparse_hybrid_200','sparse_cudss_200'):
            d=read(path/tag/'discovery.json')
            if d:
                if d['iterations_completed']!=200:
                    raise ValueError('Iteration-matched claim requires all 200 iterations')
                final=read(path/tag/'history.json')[-1]
                process=process_metrics(remote/(label.lower()+'_'+tag+'_time.txt'))
                rows.append({'backend':tag,'seconds':by_job[label.lower()+'_'+tag]['seconds'],
                    'proposal_function_seconds':d['seconds'],'process_metrics':process,'setup_seconds':d['setup_seconds'],
                    'kernel_seconds':d['kernel_seconds'],'QR_seconds':d['ideal_projection']['seconds'],
                    'ideal_coordinate_count':d['ideal_projection']['ideal_coordinate_count'],
                    'ideal_rank':d['ideal_projection']['numerical_ideal_rank'],
                    'final_b_Ha':final['b'],'final_unverified_width_mHa':final['unverified_width_mHa'],
                    'normal_map_relative_error':d['normal_map_relative_error']})
        if rows:
            reference=rows[0]
            for row in rows:
                row['b_difference_from_CPU_reference_Ha']=abs(row['final_b_Ha']-reference['final_b_Ha'])
                row['width_difference_from_CPU_reference_mHa']=abs(row['final_unverified_width_mHa']-reference['final_unverified_width_mHa'])
            comparisons.append({'case':label,'runs':rows})
    verifiers={name:read(remote/(name+'.json')) for name in ('h8_stdlib','h8_flint','h8_flint_gram')}
    if any(not value or not value['exact_endpoints_equal_to_stdlib'] for value in verifiers.values()):
        raise ValueError('All H8 arithmetic comparisons must match exact endpoints')
    sparse={label:read(remote/(label+'_sparse_corrected.json')) for label in ('h8','h10')}
    tensor=read(remote/'h10_tensor_retry.json')
    accelerated=read(remote/'accelerated_jobs.json',[])
    probes=read(remote/'probe_jobs.json',[])
    candidates=[read(p) for p in OUT.glob('*.json')]
    receipts=[row for row in candidates if isinstance(row,dict) and row.get('instance_id')==owned['id']]
    # These orchestration calls are sequential and include installation and
    # waiting for their child processes. Inner stage times must not be added a
    # second time to this total.
    measured_cloud_orchestration=sum(x.get('wall_seconds',0.) for x in receipts)
    local=[read(CASES/'runs'/name) for name in ('h10_gpu_independent_replay.json',
        'h10_library_independent_replay.json')]
    if any(not row or row['status']!='passed' for row in local):
        raise ValueError('Both local standard-library replay process receipts are required')
    prior_components=[read(CASES/'runs'/f'{n}.json') for n in (
        'h10_size_generate','h10_rotated_construct','h10_rotated_state','h10_rotated_verify',
        'h10_rotated_lower_initialize',
        'h10_rotated_lower_nonsinglet','h10_rotated_lower_prepare','h10_correlated_prepare')]
    preparation=sum(x['wall_seconds'] for x in prior_components)
    full=[]
    for tag in ('hybrid_adaptive','sparse_cudss_adaptive'):
        d=read(h10/tag/'discovery.json')
        if not d:continue
        exact=read(h10/tag/'exact_compiled_retry/complete_replay.json') or read(h10/tag/'exact_compiled/complete_replay.json')
        if tag=='hybrid_adaptive':
            lower=read(h10/'interval.json')
            proof_seconds=process_metrics(remote/'h10_replay_time.txt').get('wall_seconds',lower['lower_replay_seconds'])
            proof_kind='original stdlib lower; the prior independently checked upper is charged in preparation'
        else:
            proof_seconds=process_metrics(remote/'h10_compiled_retry_time.txt').get('wall_seconds',exact['complete_seconds']) if exact else None
            proof_kind='complete rotated upper and FLINT rational lower; independent stdlib replay additionally charged'
        process=process_metrics(remote/('h10_hybrid_time.txt' if tag=='hybrid_adaptive' else 'h10_sparse_cudss_adaptive_time.txt'))
        proposal_seconds=process.get('wall_seconds',d['seconds'])
        full.append({'proposal':tag,'proposal_seconds':proposal_seconds,'proposal_function_seconds':d['seconds'],
            'proposal_process_metrics':process,'setup_seconds':d['setup_seconds'],
            'iterations':d['iterations_completed'],'proof_seconds':proof_seconds,'proof_scope':proof_kind,
            'shared_local_preparation_component_seconds':preparation,
            'mixed_host_component_sum_seconds':preparation+proposal_seconds+proof_seconds if proof_seconds else None,
            'not_a_single_fresh_timed_end_to_end_run':True})
    summary={'best_H10':best,'all_H10_certificates':intervals,'fixed_iteration_comparisons':comparisons,
        'H8_verifier_comparison':verifiers,'sparse_system_comparisons':sparse,'tensor_comparison':tensor,
        'full_prepared_runs':full,'cloud':cloud,'remote_orchestration':receipts,
        'remote_orchestration_wall_seconds':measured_cloud_orchestration,
        'probe_processes':probes,'accelerated_processes':accelerated,
        'local_H10_replays':local,'charged_local_preparation':prior_components,
        'claims':{'H10_accuracy_met':best['target_met'],'general_scaling_proved':False,
                  'whole_DMRG_cuTensorNet_speedup_measured':False,'all_open_research_questions_solved':False}}
    with (OUT/'summary.json').open('x') as stream:json.dump(summary,stream,indent=2)
    control_lines=[]
    for case in comparisons:
        for row in case['runs']:
            peak=row['process_metrics'].get('peak_host_RSS_bytes',0)/1e9
            control_lines.append(f"| {case['case']} | {row['backend']} | {row['seconds']:.3f} | {row['setup_seconds']:.3f} | {row['QR_seconds']:.3f} | {row['kernel_seconds']:.3f} | {peak:.3f} |")
    replay_lines=[f"| {name} | {value['lower_seconds']:.3f} | {value['complete_replay_seconds']:.3f} |"
        for name,value in verifiers.items() if value]
    full_lines=[f"| {x['proposal']} | {x['proposal_seconds']:.3f} | {x['proof_seconds']:.3f} | {x['mixed_host_component_sum_seconds']:.3f} |" for x in full if x['proof_seconds']]
    report=f'''# Library integration and the H10 accuracy milestone

The original H10 rational Hamiltonian now has a locally rechecked complete
fixed-N interval of **{best['width_mHa']:.9f} mHa**. Both the orbital-rotation/MPS
upper and the full singlet-plus-nonsinglet lower were replayed with standard
Python integer/rational arithmetic. No full fixed-N determinant space enters
this construction or acceptance. The earlier 0.479085203 mHa GPU-discovered
certificate and all failed antecedent attempts remain preserved.

This is an adaptive follow-up to the frozen transfer campaign. The local-basis
MPS, transported one-/two-/three-body moments and changed numerical schedule
are explicit dependencies. It is not a success of the original timed-out
canonical random-state H10 rule, and it is not a new smaller SOS family.

## Measured library effects

SuiteSparseQR replaces a dense QR factorization of ideal-multiplier coefficient
columns. It preserves the numerical span checks, energy normalization and the
original exact acceptance requirements. This factors a constraint-coordinate
matrix, not a many-electron Hamiltonian matrix. cuDSS replaces only the shifted
normal preconditioner; the unshifted constrained solve is unchanged.

The following runs use the same host, FP64 precision, prepared arrays, zero
Gram start, 200 iterations and fixed mu=2. The setup costs are included. These
are fixed-work comparisons, not times to achieve a given certified accuracy.

| Case | Backend | Process seconds | Setup seconds | QR seconds | Iteration-kernel seconds | Peak host GB |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(control_lines)}

The machine-readable receipt includes final objective/width differences,
numerical ranks and normal-map consistency errors for assessing equivalence.
Neither an optimizer status nor a floating-point bound prediction accepts an
energy certificate.

cuDSS and CHOLMOD were also compared with SciPy SuperLU on the actual H8/H10
normal matrices. Repeated-solve timings include RHS/output transfers. The first
probe accidentally charged CUDA context initialization to its first CPU setup;
that original result is preserved and the corrected probe records initialization
separately. Its repeated-solve measurements were unaffected by that setup issue.

FLINT rational arithmetic was tested separately from compiled integer Gram
products. Every comparison below requires exactly equal upper and lower
endpoints to the original standard-library receipt.

| H8 arithmetic | Lower replay seconds | Complete upper/lower replay seconds |
|---|---:|---:|
{chr(10).join(replay_lines)}

These are paired measurements, not a broad timing distribution. FLINT rational
arithmetic is selected for the integrated experiment. The Gram-only variant
did not establish a useful improvement. Its focused tests include arbitrary
large integers, empty factors, invalid words, mixed charges and noninteger
factor rejection. The original accepting source files are unchanged.

The first integrated FLINT replay failed because the upper check had warmed a
CAR cache with Python Fraction objects, then the lower used FLINT rationals.
The retry loads the lower dependencies before selecting arithmetic and clears
that cache. A focused regression reproduces the warm-cache transition and checks
identical exact spin polynomials. The failed process and its cost are preserved.
The public replay caller also enforces the even-electron assumption required by
the integer-spin argument; all accepted fixtures satisfy it.

cuTensorNet was measured on an actual H10 MPS transfer with tensor shapes
(32,32), (32,2,63), (32,2,63), (2,2). It did not beat the existing sparse Numba
CPU transfer: about 42.7 microseconds for the CPU implementation, 120.2
microseconds for cached resident cuTensorNet, and 399.9 microseconds including
transfers. All compared values agreed to about 1e-16 relative error. This is
one contraction test, not a complete DMRG benchmark or a claim about larger
tensor networks. Its initial missing-CVXPY import failure is retained; adding
that existing transitive dependency enabled the retry.

## Complete-certificate runs and dependency cost

| Proposal | Proposal seconds | Proof seconds | Charged mixed-host component sum |
|---|---:|---:|---:|
{chr(10).join(full_lines)}

The shared local preparation contribution is **{preparation:.3f} seconds**:
integrals, exact orbital rotation, fresh local MPS, upper acceptance, nonsinglet
construction, coefficient-map preparation and transported moments. The first
proof column is the standard-library lower replay; the second rechecks the
upper and uses FLINT rationals for the lower. The latter sum conservatively
charges upper verification again. The sums combine local Mac preparation and
remote A100-host computation; they are not measurements of one fresh complete
run on a single host. Earlier failed searches, benchmarks, installation and
independent replay are additional campaign costs, not omitted dependencies.

The fresh integrated H10 proposal starts again from zero Gram/ideal coordinates;
it uses no prior winning checkpoint or FCI vector. Exact local replay is the
authority for the reported final interval. Numerical FCI was used separately
as a reference in the parent campaign, not as a construction input.

## Cost, preservation and limits

The owned A100 instance is terminated. Launch through confirmed absence was
{lifetime/60:.2f} minutes at the observed $1.99/hour rate, for an estimated
**${cloud['estimated_USD_including_boot_and_idle']:.2f}**, including boot and idle
time. The provider invoice is not available. Pre-existing instances and the
separate earlier GPU task were left untouched.

All remote job outcomes, setup failures, transfer receipts, source hashes,
installed versions and Linux process-time/memory records are retained. Linux
maximum resident sizes in the time files are KiB; local budget receipts use
bytes. Inner process durations are not added again to their enclosing remote
wall times. Earlier local research and failed H10 starts remain in the parent
all-attempt ledger.

This establishes a new H10 model certificate and working numerical-library
integrations. Broad transfer, favorable large-system scaling, competitive
complete cost and experimentally useful predictions still need evidence.

- [Exact local H10 result]({best['path']}/interval.json)
- [Machine-readable measurements]({OUT}/summary.json)
- [Additional chemistry/library shortlist]({ROOT}/research/nvidia_followup_20260915/LIBRARIES.md)
- [Broader transfer and physical-model results]({ROOT}/research/transfer_followup_20260915/RESULTS.md)
'''
    with (ROOT/'research/nvidia_followup_20260915/REPORT.md').open('x') as stream:stream.write(report)
    print(json.dumps({'H10_width_mHa':best['width_mHa'],'cloud_estimate_USD':cloud['estimated_USD_including_boot_and_idle']}))


if __name__=='__main__':
    run()
