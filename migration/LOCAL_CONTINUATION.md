# Latest checkpoint — 2026-09-13 05:38 UTC — goal ACTIVE

This goal turn made PROGRESS. Do not mark the goal complete: a general
chemistry-relevant condition for scalable accurate discovery remains open.
Keep the original goal intact. The user authorizes two warm A10s ($2.58/h
combined), CLI/SSH, and subagents. Do not tear down between batches. Never
contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

## Live work: H10 only

Instance A, 129.159.32.200, is idle after completed H8 source and proof runs.
Instance B, 146.235.200.232, is running H10 SU2 SCS600. Exec session55605,
Python PID16627, remote output:
/home/ubuntu/spectra-spin-discovery-v2/results/spin_runs/h10_scs.
Construction took681.690s, with10,040,532 map nonzeros. Solver cap remains600s.
The original whole-process timeout1500s was too short for build+solve+export;
root installed a bounded export guard with a2100s whole-process ceiling.
Guard exec session28736, remote script spin_export_guard.py, log
/home/ubuntu/spectra-spin-discovery-v2/results/spin_runs/h10_export_guard.log.
The guard suspends ONLY timeout PID16626, supervises the Python child, and
resumes the timer when the child completes or reaches2100s. Do not manually
resume or kill that timeout during export. Local guard copy:
results/lambda_runs/scs_refinement/h10_export_guard.py.

Running cloud code is frozen v2:
SHA7b45026d6fe407fa00e9cef636fa317e1bc9a8658d249654d26d4ec5057fcf10.
Root has a newer LOCAL v3 source, NOT deployed over the running code:
results/lambda_runs/scs_refinement/spin_discovery_v3_source.tar.gz,
SHA8058408841be83b73ad9383bbc9fe4f155de8be33723152e40aa0bff58f063ac.
No agents remain assigned active work. All completed; root owns their files.

## Exact SU2 discovery backend

Read research/certificate_scaling/spin_irrep_frontier.md and campaign
results/certificate_scaling/spin_irrep/campaign.json.
Modules spin_basis.py, spin_parity.py, spin_invariant_discovery.py implement
complete rational highest-weight decomposition of the actual chosen word
families, exact compatible parity intersection, and all cross-copy Gram terms.
Local spin pattern dimension is at most8. Every descendant/rank gate is exact.
Final certificate retains ORIGINAL H; rounded-H symmetry defects are charged
in its full residual. Coefficient-l1 residual objective is not SU2 invariant;
do not claim equal finite-residual optima to the unrestricted formulation.

Invariant Gram entries H4/H6/H8/H10 are5088/56996/320192/1224400, compared
with17776/217268/1206464/4674900 previously. The corrected basis work counts
are in decomposition_ladder_corrected_counts.json. Earlier counts missed
validation actions, but dimensions and whole-run times were correct.

Independent intervals:
- H4 SCS30:7.898699e-7 Ha PASS; source6.298s, replay.611s.
- H6 SCS180:6.81108637e-5 Ha PASS; source209.352s, replay6.906s.
- H6 Clarabel120:7.189853963 Ha FAIL.
- H4 Clarabel60: coefficient width.001694 FAIL; spectral width.000375914 PASS.
- H8 SCS300: source555.279s,2500iterations, coefficient width6.276704644 FAIL.
  Build147.973s, solve304.662s, export102.644s. Negative eigenvalue mass.087609.
  Spectral residual proof improves width to.699480442455, still FAIL.
  Proof90.081s plus independent spectral interval replay93.802s.
  Local h8_scs/ and h8_scs_spectral/ are complete, hash-matched to remote.
  CertificateSHA34ead5365879bcb7f2ec02219efb8124235c7b52ecabfaa01e199ef49f663aad.
  ProofSHAe87f83f213be28f3254948d54b86b56b4f43ae46d8c3431c650727f541327793.

The older unrestricted chosen mixed-cubic H4/H6/H8/H10 ladder still all passes
.0016 Ha; do not replace it with these failed finite-budget SU2 trials.
That dictionary has all mixed cubics but selected local pure triples only.
FCI upper discovery remains an exponential validation cost.

## Adopted exact map optimization

Profiling exposed repeated canonical/product/Hermiticity/Fraction operations.
Agent's memoization suggestion was redundant: word_product already has LRU.
Its first fused implementation was9% slower and was rejected. Root rewrote
spin_gram_columns.py with precomputed channel adjoints and exact canonical
monomial adjoint permutation signs, retaining exact Hermiticity checks.
All28,894 H6 columns, indices and word-pair counts match the independent CAR
reference. Cold-cache times16.024s reference versus9.452s revised (one local
comparison, not a total solver speed claim). H4 exhaustive and rational
cross-copy tests also pass. Backend now imports this helper for FUTURE runs.
Column benchmark: results/certificate_scaling/spin_irrep/column_builder_benchmark.json.
52 focused tests passed after wiring; focused_tests.log. No full-repo regression.

## Interacting parent condition beyond hidden quadratic controls

Read projected_product_parent.py, test_projected_product_parent.py and .md.
Root replaced agent drafts with incorrect density signs and mirrored tests.
Public recognize infers the full connected hopping path and a real sign gauge,
transforms EVERY original term by that signed mode permutation, and checks
its inverse exactly. Recurrence alpha0=A0, beta=t²/alpha,
nextalpha=A_next-beta recovers positive rational local parent weights.
Require density/boundary equations and actual CAR reconstruction of the
ENTIRE H as c+sum(g B†B). A common fixed-N projected-product null state proves
E0=c. Normalization uses elementary symmetric DP, no occupation enumeration.
The inverse signed permutation transports its fermionic wedge signs.

Six focused tests include independent N2 Fock action and wedge permutation
signs, rational weights, all sectors' recognition, mutations, CLI under-S,
and a branching-star counterexample. The unchanged parent factors on a
four-mode star have N2 determinant4 and satisfy(H-2)^3=3(H-2), with minimum
2-sqrt3>0. Thus treewidth1 alone does not extend THIS common-null ansatz.
This is not a theorem that all certificate methods fail on trees.

Ladder M4..256, N=M/2: exact width0. The current inferred-path checker takes
.1457s atM256,255factors,1020factor monomials,32768DP updates,508 maxDPbits.
Older natural-order core timing.0913s is separately preserved. Current results:
results/certificate_scaling/projected_product_parent/inferred_path/ladder.json.
All four molecular fixtures REFUSE. The state is non-Gaussian (N2 Pluecker
violation), but the parent mechanism is related to established RK/SMF work;
no general-physics novelty or molecular applicability claim.

## Next actions

1. Wait for H10 source and guard terminal status. Download, hash-check, and
   independently replay original-H interval. A cap/failure is not a pass.
2. Use completed H8/H10 diagnostics to decide conditioning/cone work; don't
   assume fewer Gram entries means better accuracy or total runtime.
3. The faster map helper is frozen locally for subsequent experiments.
4. Keep both authorized hosts warm and goal active. Update this checkpoint
   after new outcomes. Do not repeat completed jobs or inflate claims.

---

# Latest checkpoint — 2026-09-13 04:57 UTC — goal ACTIVE

Previous goal turn is PROGRESS: H10 crossed the accuracy threshold under
independent exact replay; new conditioning experiments and exact spin controls
completed. DO NOT mark goal complete: no general chemistry-relevant structural
accuracy/discovery theorem. Keep full objective intact.

## Authoritative current state

ALL remote and local research jobs terminal. At04:56UTC SSH ps on both hosts
found no spectra-venv Python process. A129.159.32.200 and B146.235.200.232 stay
warm by user authorization ($2.58/h combined). Do not teardown or duplicate
old jobs. Never contact archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

H10 accelerated cold1800 completed2475iterations. Independent local python-S
interval is .001157618253846804 Ha PASS (lower-12.369457074577289,
upper-12.368299456323442). The H4/H6/H8/H10 frozen accuracy ladder now passes:
1.530023922e-7 / 5.499439284e-5 / .000750585102643 / .001157618253847.
Scoreboard results/certificate_scaling/active_space_ladder_latest.json.
Full chosen mixed-cubic dictionary Gram entries17776/217268/1206464/4674900;
export factors499/1116/2574/4810. This family has all mixed cubic words but
only selected local pure triples, not literally every possible cubic word.
Finite passing ladder is NOT an asymptotic accurate sparse-discovery theorem.

H10 artifacts:
- results/lambda_runs/scs_refinement/downloaded/h10_evd_1800/
- results/certificate_scaling/scs_refinement/h10_evd_1800_spectral/
- results/certificate_scaling/scs_refinement/intervals/h10_evd_1800.json
Source2222.9477s +spectral255.1127s +independent182.8796s. OriginalFCI
upper-witness discovery still costs exponentially and is only validation.
CertificateSHA397eb5bbed0dc09b24169d08307f448aea1d47749828d7d32f1470a912cdf7c0.
SpectralSHA24756bfae23b26704bda8c8cf83ab1e2dc38e57c508b2ff9948a2ee30354bc58.
Twenty new output files hash-match remote originals; whitening_download_hashes.json.

## Completed conditioning experiment

Read research/certificate_scaling/polynomial_basis_conditioning.md.
Root modified polynomial_gram_contraction.py and commutator_dictionary.py:
optional --basis-condition whiten (requires --map-backend contraction).
Numerical C'=C*W.T; exporter computes sqrt(Q)*W then rounds effective factor
rows in ORIGINAL exact polynomial basis. Original CAR residual pays all error.
No automatic rank truncation; numerical rank deficiency refuses.
Whiten separately in connected components sharing EXACT monomial support to
preserve structural zeros. Cross-component physical Gram variables remain.
Whole-block SVD produced numerical fill: H4map2.3927m instead167104; archive
retains transforms/receipt/cert but transient whole-block code was superseded.
ComponentSVD H4map167104, build.604s/solve2.393s/total3.235s, exactwidth1.23319e-6.
Depth0whitenedH4control width.003566801899 still obeys prior exact obstruction.

H6 component whitening has43740Gram,9419rows,4839646mapnnz:
Clarabel120 exactwidth .0560815951320, build11.639/solve137.527/total161.712s;
SCS180 exactwidth .210586928014, build11.535/solve231.223/total254.287s.
BOTH FAIL. Better than unwhitened finite runs but no cone optimum proof.
No further H6 Krylov jobs are running. Avoid unmotivated longer repeats.
Output dirs commutator_krylov/{h4_whitened_clarabel,h4_component_whitened,
h4_depth0_whitened_control,h6_component_whitened,h6_component_whitened_scs}.
Frozen cloudsource /home/ubuntu/spectra-whitening onA, source archive
whitening_source.tar.gz SHA2ad69feb4bb1ace7496db5fd337c65a2578374810bf7aabacbf52eb69573057b.
Root rejected agent operator_pricing's false claims: diagonal commutator
CAN add directions (different weights per word), H6 basis not shown worse
condition than H4, and physical dual=D*y_scaled is exactly correct. Exact
new test captures diagonal newdirection. Do not repeat those false claims.

## Next structural attack: SU(2) operator multiplets

Read research/certificate_scaling/spin_irrep_frontier.md. Backend NOT built.
New tested algebra module spin_multiplets.py: highest-weight checks, divided
power descendants, exact raising/lowering/Sz identities, invariant squares
and ORIGINAL-format integer factor export. Integer repetition weights:
j0[1], j1/2[1,1], j1[2,1,2], j3/2[3,1,1,3]. Correspond to common positive
multiple of inverse-binomial weights. Tests include independent occupation
matrix commutators and wrong-weight rejection. Creationhalfspin square is
2-nup-ndown, NOT nup+ndown (initial agent example was corrected byroot).

Original H fixtures have tiny exact spin defects from rounding. Used EXISTING
experiments/marginal_spin_reduction.symmetrize/commutator/spin_operators.
Exact Casimir projection Hs and [S+,Hs]=[Sz,Hs]=0 verified in python-S.
||H-Hs|| <= coefficient L1: H4 2e-12; H6 5.4e-11; H8 2.36e-10; H10 4.74e-10.
Artifacts results/certificate_scaling/spin_irrep/{hamiltonian_symmetry.json,
projection.json,h4,h6,h8,h10}. Projected fixtures are proposals ONLY; final
certificate must contain ORIGINAL H and pay complete residual. Do not silently
replace benchmark H or upper witness with Hs.

Suggested next bounded implementation:
1. Build highest-weight multiplicity bases of actual dictionary families.
   Spin ladder action preserves each spatial orbital's separate creation/
   annihilation counts, so degree<=3 kernels can be grouped into small local
   support patterns. No full Fock space. Validate completeness/dimensions.
2. Intersect GF2 parity-generator span with masks equal on each alpha/beta
   orbital pair. Cannot just filter individual generators: commuting linear
   combinations may otherwise be lost. Preserve number charge and spatial sectors.
3. Use invariant descendant sums for each multiplicity Gram matrix. Preserve
   cross terms between equivalent copies. Measure map fill and total discovery
   cost; fewer PSD blocks alone may be slower. Do H4 exact original-H replay
   before bounded cloud H6 etc. No generic fixed-degree accuracy guarantee.
4. Residual-l1 objective is NOT SU2 invariant; averaging is not proof of
   equality of finite-residual optimization objectives. Preserve exact gates.

Existing dictionaries are experiments/marginal_coefficient.dictionaries;
parities experiments/marginal_molecule_stress.parity_generators. Baseline
adaptive uses alpha-charge partition, not SU2. Number-ideal scalar projection
is a further choice; no new backend or highest-weight nullspace exists yet.

41focused tests pass in .852s (whitening_and_spin_tests.log); separate4spin
controls pass python-S. No full regression claim. All agents completed, root
owns modifications. Current fullsource spin_frontier_source.tar.gz SHA
c9d78567b07864f1061eb0f2cb9d5f014ead3d2b9505df3945dd0f749536dcbf (347files).
Campaign index updated no live jobs; SCS_REFINEMENT_RESULTS.md starts with
latest passing milestone. Both hosts remain warm for next authorized work.

--- Older checkpoints below are historical where superseded ---

# Latest checkpoint — 2026-09-13 04:31 UTC — goal ACTIVE

Continue the existing goal; DO NOT create a replacement or mark it complete.
No chemistry-relevant general sufficient discovery condition has been proved.
Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.
Both A10s remain warm by user authorization ($2.58/hour combined). Use SSH.

## Immediate live work

ONLY B H10 accelerated cold1800 is live:
- Host146.235.200.232, PID13348, exec session92363.
- Source /home/ubuntu/spectra-evd-refinement, immutable evd_source.tar.gz.
- Output /home/ubuntu/spectra-evd-refinement/results/h10_evd_1800.
- Log /home/ubuntu/spectra-evd-refinement/h10_evd_1800.log.
- Started ~04:02UTC; at04:27 runtime25min, CPU100%, RSS4.4GB.
- SCS1800s, overall2300s, external2700s; exact export follows solver.
- EVD preload only this process; same verified shim/library as below.
- Do not duplicate. After completion download, hash-check, run spectral
  residual proposal on remote B then independent local python-S interval replay.
  Spectral run is another ~260s; all costs count. Reference H10 path remains
  results/certificate_scaling/cubic_precision/reference_h10_full_checked/h10/upper.json.
A129.159.32.200 is idle and intentionally warm. All H6 jobs terminal.
All local replay/test/sweep handles are terminal; avoid polling obsolete sessions.

## New verified molecular results and failures

H10 accelerated600 independent spectral width .00386819649650069 Ha (FAIL),
improved from stockcontinuous1200 .007893092560296. Lower-12.372167652819943,
upper-12.368299456323442. Total source1035.496s plus spectral~258s plus
independent replay. 825SCSiterations. Full cubic discovery still4674900Gram
entries;4907exportfactors,2656820factor nonzeros. Numerical convergence win,
NOTstructuralcompression. See SCS_REFINEMENT_RESULTS.md and campaign_index.json.
Source certificateSHA e73a14d4cad8fe48f03a0737c697c9f5533fd44cbfd97dad7f8233c8a8a7b1aa.

H6 depth1 Krylov all FAIL (exact independent interval widths):
fullh1 Clarabel120 6.63125026986; diagonal Clarabel120 .353743685967;
diagonal numerical contraction SCS180 18.1236494791; same rowconditioned
SCS180 4.24922911422. These are solver outcomes, NOT cone impossibility proofs.
Directories results/certificate_scaling/commutator_krylov/h6_{one_step,
diagonal_one_step,diagonal_scs180,diagonal_scs180_conditioned}.

Root polynomial_gram_contraction.py: exact monomial CAR map then floating
batched contraction, no threshold clipping, exact exporter unchanged. H6
construction153.612->14.516s, wordpairs3190092->200340, same43740Gram and
4839648mapnnz. Whole H4 3365964 map coefficients compared, error<=1.776e-15.
H4 contraction Clarabel width6.51119e-5, build.681s,total4.466s. Conditioned
SCS30 H4 width.000106215166. Basis coefficient SVD diagnostics fullnumerical
rank: H4condition124-416,H6condition145-350; not exact rank/convergence proof.
commutator_dictionary.py new --map-backend exact|contraction, --solver,
--row-condition; physical saved dual rescales correctly. No orthogonalization
implemented. Next numerical move should diagnose solver formulation before
more unmotivated long runs. The older exact H4 dual only excludes depth0.

23 new downloaded artifact files hash-match remote originals; receipt
results/lambda_runs/scs_refinement/new_download_hashes.json. Logs separate.
Latest frozen source structural_source.tar.gz SHA
782a6ae3a73c82fd5e47176e6cdc6c3015533b9c9feb103325d130fa5561eee0.
H6 conditioned source conditioned_contraction_source.tar.gz SHA
f8252889ae9f191b05671e3041cdf67d88f2be4efdfd817e547f4fc96ac3adda.

## New structural control: matching-CZ-dressed quadratic recognition

Read research/certificate_scaling/dressed_quadratic_structure.md.
Root completed/corrected agent drafts in dressed_quadratic_recognition.py,
dressed_quadratic_certificate.py and both tests. Agents are done; root owns.
- Recognize H=UqUdagger from H only for an unknown matching CZ transformation.
- Exact canonical rational/Hermitian/number gates. Graph remains connected
  after each vertex deletion; BFS solves each dressing column up to a bit.
- ALL-PAIR symmetry constraints determine at most two complement candidates.
  Check matching and regenerate entire H exactly. External control sets only
  are complement invariant. Pure quadratic validated shortcut supported.
- Compact factors: positive spectral part uses annihilation p=a =>a†a;
  negative uses creation p=a†=>aa†. Exact rounded hole-row norm sets qtrace.
  b=c0+mu*N-qtrace; number_multiplier is mu IDENTITY (not number!).
- Full H unchanged, factors dressed to cubic, integer original verifier export.
- Independent upper is rational occupied C with gamma=C(CtC)^-1Ct, transported
  by verified U. No FCI or sector vector in discovery/acceptance.
- CLI exact replay python -S -m research.certificate_scaling.dressed_quadratic_certificate
  --witness .../witness.json --out .../independent_replay.json.

Permuted two-hub family sweep M4,8,16,32,64 at half filling all exact-S pass.
Widths4.8951e-9,5.5913e-8,3.4718e-7,1.6384e-6,7.0554e-6. M64 has64factors,
293111-byte whole witness,9.933s proposal+embedded replay+extra replay;
separate-S replay5.220s. Artifacts results/certificate_scaling/dressed_quadratic/.
Both witnesses O(M²) coefficients, bounded-degree SOS replay O(M³) arithmetic;
precision/bit growth count. Fixed floating proposer no universal accuracy or
termination guarantee. Family is hidden-free/integrable, NOT molecular success.
Actual H4/6/8/10 fixtures all explicitly refuse unsupported quartic terms.
Molecular transfer report saved. Do not mark the main goal complete for this.

Independent agent independent_audit read final files and ran10tests; no concrete
soundness bug found. Root focused campaign32tests pass in .781s; separate
recognizer6tests pass python-S. Tests now genuinely include independent full
M4 occupation-mask CZ identity, fixedN0..4 numerical spectra, malformed and
rank-deficient refusals, and exact residual payment. Earlier agent compile-only
claims/incorrect test drafts were rejected and replaced by root. No full repo
regression claim. Literature precedents Fendley1901.08078 and Jones/Linden
2107.02184; our note makes no novelty claim about hidden-free physics.

--- Older checkpoints below are historical where superseded above ---

# Latest checkpoint: exact commutator obstruction, H10 refinement, LIVE experiments

Goal remains ACTIVE. Do not mark complete. Read research/certificate_scaling/
SCS_REFINEMENT_RESULTS.md and commutator_krylov.md. Index results/certificate_scaling/
scs_refinement/campaign_index.json. Both A10s warm by user authorization ($2.58/h total).
Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

NEW completed: continuous1200s H10 exact spectral interval .0078930925603 Ha
improves previous .02334631, still FAIL. Restart2x600 gives .20824579599 FAIL.
Both root python-S independently replayed. Source and all18 downloaded file
hashes matched. Continuous total1638.08s +spectral262.81s +independentreplay;
restart1734.72s +spectral279.64s +independentreplay. Rawgap NOTcertificatequality.

Exact H4 creator-channel dual obstruction PROVED by root stdlib checker
commutator_dual_witness.py: y(H)=-3.6705580398652247, physicalgroundlower
-3.667000108966591573 =>restrictedL1 coneerror>=.00355793089863313>.0016.
Fullcharge64,64,28,28,64 exactPSD,443body2idealconstraints,|y|<=1. 85,170bytes.
Allnumericalgeneratorscontainedinrebuiltcone exactchecked. Refuses Krylovsteps>0.
Exactdual artifacts commutator_dictionary/h4_exact_dual. Does NOTexcludeother
residualbounds/dictionaries. Rootcorrected agent's offdiagonal factor2 error.

NEW diagonal_driver option: [diag(h1),p] exactsameword energyweight; H4
exactwidth6.51245170e-5 PASS,8304Gram,167104map,6.605s total. OriginalfullH1
uses274580map/10.335s. No H6diagonalrunyet; nextboundedcandidate oncecurrent
H6full-driver finishes. Frozen current H6 source doesnotinclude diagonaloption.

Root implemented one_body_steps optional H1 commutator Krylov extension:
H4_true_one_step exactwidth4.89771659e-5 PASS. Gram8304 vsfull17776 but
map274580 vsfull10256; 10.335s vsfull~7.26, NOTperformancewin. H-only, no
reference/dualfactorsinput. Agenth4_one_step_v2 falselynamed oldrun; annotated.

20focusedtests pass (checkpoint3, adaptive7, commutator6, dual4).
scs_checkpoint.py saves actual canonicalx/y/s for inaccurateiterates, explicit
identity/load/refusal. Records declared backend override libraryhashes when set.
H4resume exactwidth3.130e-6; same/nooverride older checkpont identity unchanged.

PSDkernel profile95%conetime. SCSbundledOpenBLAS0.3.15Prescottsinglethread
actual8x725 dsyev2.275s vsdsyevd1.151s. RootcontrolledSciPySkylakeXone-thread
8x725evd.4488s, ev1.7593s. Rootreplaced baduncontrolledagentbenchmark; original
receipt markedEXCLUDED. Allkerneldata scs_refinement/{projection_gate,alternate_projection_gate}.

Rootcompleted scs_evd_preload.c: isolatedLP64dsyev->scipy_dsyevd via lazy
RTLD_LOCAL|DEEPBIND, bothworkspacesallocated. Initialagentworkspaceoverflow
anddirectlinkOpenBLASsymbolcollision debugged; failedlogsretained. Stock/shim
4x4+20x20controls125iters objectiveagreement<3e-14,504interceptcalls/noerrors.
H4actualmolecular-Sreplay width3.8919e-7,134442calls/noerrors. Nevergloballypreload.
ShimB /tmp/spectra-evd-gate/libspectra_evd.so SHA8e61ca79a56199d2766419136425a9298800dba22192386f4adbf2d6fe3b6bcb.
SPECTRA_EVD_LIBRARY=/home/ubuntu/spectra-venv/lib/python3.10/site-packages/scipy.libs/libscipy_openblas-c128ec02.so

LIVE JOBS (poll SAME sessions, do not duplicate or teardown):
* B146.235.200.232 H10cold600s fasterbackend, exec47951, PID12589,
  /home/ubuntu/spectra-evd-refinement/results/h10_evd_600,
  log /home/ubuntu/spectra-evd-refinement/h10_evd_600.log.
  600solverseconds,1200overall,external1500timeout;checkpointandrawsaved.
  frozen evd_source.tar.gz SHA7aa4db0a7126cd8fcae3fdbbba2c25005cc8dd63ba74dafe2eaf5f0dbdfee1cc.
* A129.159.32.200 H6creatorKrylovdepth1, exec21850, PID11412,
  /home/ubuntu/spectra-krylov-refinement/results/h6_one_step,
  log /home/ubuntu/spectra-krylov-refinement/h6_one_step.log.
  120sClarabelsolve,external900timeout. frozen krylov_source.tar.gz
  SHAd5685f03ad27cd2652a02844770e71984b30b9fe693dbe15afe7fafd27cd81e9.
Need download/results/exactintervalreplay aftercompletion, accountallcosts.
IfH10needsstrongerresidual use wedge_spectral_bound remotely then root-S
cubic_interval_replay. ReferenceH10path in previous report/index unchanged.
AlloldSCScontinuous/restart/projection/localreplayjobs terminal. Agents stopped;
rootownsallcurrentmodifiedfiles. Do nottrustagentclaimedexperimentswithoutreceipts.

--- Prior checkpoint superseded where above differs ---

# Latest checkpoint: cubic precision campaign completed; goal ACTIVE

Authoritative report: research/certificate_scaling/CUBIC_PRECISION_RESULTS.md
Machine index: results/certificate_scaling/cubic_precision/campaign_index.json
Do not mark goal complete. No general accurate structural discovery condition.
All campaign/local replay processes terminal. Both A10s intentionally warm.
Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

NEW best uncompressed exact intervals (python -S replay, same frozen rational H):
* H4 width1.530023922e-7, original coefficient bound, Clarabel local7.259s.
* H6 width5.499439284e-5, SCSconditioned259.309s +spectral residual.
* H8 width7.505851026e-4, SCSquotientconditioned721.563s +spectral residual.
* H10 width.02334631083 FAIL, SCSquotientconditioned1215.936s +spectral.
  Only175SCSiterations/673.245solve seconds. Cannot infer cone-optimum gap.
  OriginalL1 width.03221293, Gersh.02566182. Needs convergence diagnosis.
Final exact Fraction intervals/hashes in intervals/final_h{4,6,8,10}.json.
H6originalL1already .00025153. H8Gershalready .00148134. Spectral strengthens.

Root changes: adaptive_block_discovery.py raw snapshots, precisiondiagnostics,
row/columnconditioning, number-ideal linear quotient and identity anchoring,
optionalexplicitdual. ClarabelH6primal2+dual1all360stimeout, saved logs.
H6quotientSCS worse than nonquotientconditioned; don'tclaimquotient universal.
H8unconditionedwidth.026285 thenconditioned.002974 beforestrongerreplay.
Full dictionary Gram entries H4/H6/H8/H10:17776/217268/1206464/4674900.
Full SOS rows499/1116/2574/6233; discoverystillfullandallcostcharged.

wedge_residual_bound.py exact fixed-N exterior k-body Gershgorin theorem.
ROOT REPLACED wedge_spectral_bound.py accepts integer factors, computes exact
Gram residual, canonicalphase(-1)^(k(k-1)/2), real k-tupleindices, minimum
acrosscomponentsincludingzeros, thencomb(N,k)lift andmaxoldL1/Gersh.
5tests +jointaudit pass. Oldagentbase wedge_spectral/receipt.json INVALID
(collapsedindices); use h6/h8/h10_conditioned subdirectory witnesses ONLY.
Added residual factors298/696/1350 andbytes95108/554403/2190186 counted.
Root extracted shared verified_residual() in experiments/marginal_symbolic.py;
verify publicreceipt unchanged, everyvalidationgate retained. Gram symmetric
half-dot products +contiguouscolumns. tests17symbolic/orbit+20researchpass.

Root wedge_multiplier_repair.py sparseLP keepsH/b/factorsfixed, repairsX.
H6Gershwidth .00008440851; H8unconditioned .00991820FAIL. OriginalL1may
becomemuchworse fromcancellingX; use namedwedgereplayonly. TestmutatesX,
notH, andrestoreszeroexact. Sourceinitialagentdraftwasreplacedbyroot.

Compression AFTERdiscovery: H6unconditioned238rows/33400nz width.00131997,
H6conditioned252rows/34216nz width.0003966475. H8conditioned576rows/
191936nz +696residualfactors width.0008204347; originalcert2125203bytes,
proof559994bytes. 3H8trials115.952s +final-Sreplay38.152s. All full
source solve costs still counted. Artifacts factor_structure*, all exact.

Reference uppers reconstructintegralsfromfrozenH, correctsigns,nofactor4,
correctFCIinterleavedparity. Full H8reference2468nonzerostates; H10=31752
of63504selectedpairs, remainingroundedzerosnottruncation. Streaming integer
upper O(K*terms), nocache/K²matrix; oldOracle4096gateunchanged. H10exact
upper -12.368299456323442,227026800checks~16s. Validation-onlyFCIstill
exponential cost, forbidden asdiscoveryinput. All finalintervals rootreplayed.

GPUkernelgate: A10slower8x368/16x112. For8x725GPU.22937s vsCPU.34981s
(1.52x),8x225slower. Isolated float64kernelonly,notendtoendsolverspeedup.
SCSresume nextread-onlyfinding: CVXPY1.9.2cachesONLYoptimal, NOT
optimal_inaccurate. Read scs_resume_design.md; verifyremote1.6.5source
beforeimplementation. CurrentrawNPZhasphysicalx/rem/dual/Gram, notfull
canonicalSCSx/y/s. NeedexactproblemidentityhashandpersistSCSstate for
boundedrefinement; defaultwarm_start alone won'tresume inaccurateiterates.
Secondroute note positive_amplitude_structural_route.md: existing
stoquastic_chain_certificate.py alreadyknownspinparentcontrols256sites;
fixed-Nfermion+largerlocalJastrowDP is proposedextension, notnewcompleted.

Hosts: A129.159.32.200 idfab09ca8165b41d7bfc5649f413dda3f;
B146.235.200.232 idb55288de9761445f8c589c21829f853c. SSHubuntu,keydefault.
Both$1.29/h ($2.58combined),keptwarmbyuserauthorization. No teardown.
Remotev1 /home/ubuntu/spectra-precision, v2spectra-precision-v2,
PYTHONPATHfallback /home/ubuntu/spectra-structural. Venvspectra-venv Python3.10
cvxpy1.6.5,SCS3.2.8;systempythonTorch2.7CUDA. No new cloudjob currently.
42downloadfilehashesmatched. source.tar.gz,source_v2.tar.gz,final_source.tar.gz
+final_manifest.json in results/lambda_runs/cubic_precision. Preserveoldwork.

--- Prior checkpoint (superseded where above differs) ---

# Active structural-scaling goal — current checkpoint

LATEST COMPLETED CONTINUATION: hidden-basis campaign. Authoritative new report
research/certificate_scaling/HIDDEN_BASIS_RESULTS.md. Goal remains ACTIVE;
do not call completion for positive controls. All processes terminal, both
existing A10s warm intentionally. Never contact archived tasks.

* Root implemented hidden_density_basis.py: exact CAR commutator map on
  real-symmetric one-body A, SVD/eigenbasis proposal, bounded rational
  reconstruction, exact orthogonal orbital transport via exterior-square
  tensor map. Input H only, no generating rotation/factors/state. Hidden
  known interacting ratio family M4/6/8/10/12/14/16 all exact [0,0]. M16:
  15 inner factors, 1,053,465 bytes including dense H, 22.97s discovery,
  1.11s transport replay on remote CPU. Known solvable controls, not novel
  physics/general chemistry. Exact bounded-denominator recovery may refuse.
* All canonical molecular H4/H6/H8/H10 real commutant nullities <=3 proved
  by modular rank minors, insufficient for density bases (need M8/12/16/20).
* NEW quantitative obstruction commutant_obstruction.py and
  approximate_commutant_obstruction.md. Project quartic V off all one-body
  wedge lifts L(A) (quartic part Nhat Q(A)); root exact CAR implementation
  orthogonality-checks symmetric AND antisymmetric lifts. Input must be real
  rational Hermitian; imaginary projection coefficients then vanish, so
  computed projection is full Hermitian projection. Full P is U(M)-equivariant
  and takes density-diagonal V0 to density-diagonal P(V0). Arbitrary X drops
  out via P(L(X))=0, including off-diagonal X. Joint auditor raised both
  loopholes, then withdrew after these proofs and explicit gates/checks.
* Full Hermitian commutant has real-symmetric and imaginary-antisymmetric
  Gram blocks; square-root-free metric diag(1,2). Exact integer Bareiss
  inertia at t=1/10000 gives k below t; residual Frob² >=(M-k)t/4.
  This rules out ANY complex-unitarily-rotated density-density quartic V0
  PLUS arbitrary one-body-multiplier number ideal with coefficient-L1
  residual <=.0016Ha. Quotient H4/H6/H8/H10 k=4/4/4/8 respectively; exact
  distance floors .01/.014142135623/.017320508075/.017320508075Ha.
  H10 without quotient is .02Ha; do NOT reuse that for quotient case.
  Frobenius is quartic wedge coefficient norm, NOT fixed-N operator norm or
  ground-energy error. No obstruction of general SOS or higher-body ideals.
* hidden_campaign_replay.py run python-S:7transport+4rank+4quotient proofs
  all replayed,107.64s.16focused tests pass. Results/validation:
  results/certificate_scaling/hidden_density_basis/validation.json and tests.log.
  8initial remotejobs all exit0,22resultfile hashes matched. Initial archive
  results/lambda_runs/hidden_density_basis/source.tar.gz, later proof archive
  final_source.tar.gz + final_manifest.json. Earlier H10 real-inertia run
  failed only decimal hugeint hashing; fixed hex hash, accepted complex and
  quotient results replay with final code. All jobs stopped; hosts staywarm.

NEXT ACTIVE ATTACK (supersedes old hidden-basis TODO below): noncommuting
molecular certificate precision/conditioning, then factor compression.
H6 full cubic prior numeric lower -6.3336177177 vs exactcorrelatedupper
-6.3330586262 would fit.0016 target, BUT numericlower NOTcertified, exactexport
was -6.3410288521, repair -6.339904088. Diagnose coefficient equality error,
PSD clipping and multiplier cancellation; save raw cone matrices; test
better conditioning and solver precision/alternative installed backends on
warm hosts. Need exact certificate toaccept. H8 prior numericfloor already
-9.301819 vsupper~ -9.252, so need representation/solve improvements too.
Root research/adaptive_block_discovery.py owns implementation, no source
factors/upper allowed discovery; upper only evaluation. Check available
solvers before choosing dependency. Allfullmaps/pricingcostmustcount.

--- Previous campaign checkpoint (historical next-attack text superseded) ---

The user explicitly asked to CREATE A GOAL and attack structural conditions
from all angles. An ACTIVE goal now exists (unbudgeted) to discover and
rigorously validate a nontrivial efficiently checkable structural condition
for accurate fermionic certificates with controlled discovery/verification
growth, evaluated against a growing active-space ladder. Do not create a
second goal, mark it complete for a manufactured solvable control, or reuse
historical no-active-goal statements below. Chemistry scaling remains OPEN.

Authoritative new report: research/certificate_scaling/STRUCTURAL_ATTACK.md.
The current research batch is COMPLETE; no local/remote research process live.
Two A10 hosts remain warm intentionally at$2.58/h combined. No newlaunches.
Never contact archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

NEW EVIDENCE
* 24 remotejobs:20quadratic/adaptive+4cubicladder.23exportedcerts,1H10cubic
  300s timeout no certificate.95resultfilehashes+4archivesmatched.23remote
  intervals and1fixedfactorH6idealrepair exactlyreplayed via python-S.
* H4/H6/H8/H10 genuine straightchains at1.4Angstrom,STO3G, canonicaland
  occupied/virtualblocklocalized fixtures. Bestcanonical widths .0000486761,
  .0068454618,.0269978844,.0768817055Ha. Target .0016 TOTAL Ha; onlyH4passes.
  Originalsquare/rectangleH4newfullcubic widths .0000448304/.00000799Ha,
  bytes97122/78868. No H6+ chemicalaccuracy scaling achieved.
* Fullquadraticsymmetry improvescost. H10entries233000->44300. Fixed5round
  smallPSDsupportsearch deteriorates badly onlargersystems. No representation
  impossibility claimed. Actualalgorithmrootowned adaptive_block_discovery.py
  sparsevectorizedmaps+dualpricedprincipalblocks,retainsallhistoricalsupports;
  --full runsfullchargeblocks; --symmetry exactalphaU1+GF2; notSU2newbackend.
* Everyorbitalgraphcomplete, treewidth3/5/7/9, canonical+blocklocalized.
  coefficientL1deletionbudget.0016cannotremoveanyedge; minimumedgecostfor
  localizedH10~1.435Ha. Scopegraphdefinitionandcoefficientnormonly.
* Exact locality_dual_counterexample.py disproves gap-only=>restricteddual
  pricingdecay. H=gappednearesthopvacuum, restrictsinglemonomialsquares,
  NO numberideal; exactendpointRayleigh−1 atdistance3/7/11/15. Rootreplaced
  agentfakeconcatenationCARcheck. Use locality_dual_exact.json only.
* NEW stoquastic_chain_certificate.py: positiveJastrowlocalamplitude ratios,
  local2x2PSDdecomposition+exactcyclicDPmin/max. Designednoncommuting spin
  parentfamily8..256spins all[0,0].256certificate10900bytes,discovery.535s,
  replay.027s. Notmolecular or genericfermionicresult. Numericoptimizer not
  certifiedbestparamproof. Fourindependentmatrix/DP/controltests pass.
* NEW fermionic_ratio_chain.py: H-only exactpatternrecognizer forknown
  biasedexclusion/XXZ-typeopenchain, r inferredboundarycoefficient. Hquartic,
  cubicQ_i factors, fixedNψ=r^(sum i ni); normDPpolynomialbitcost.256orbitals
  N128:255factors1020coeffs81788bytes,discovery.068s,replay~1.4s,[0,0].
  Exactsmallsectoraction/noncommutation/refusaltests. This isdeliberately
  solvableparentfamily positivecontrol; doNOTcallgeneralchemistrybreakthrough.
* 10focusedtests pass acrossadaptiveblocks,stoquastic,fermionicratio.
  Validationresults/certificate_scaling/structural_intervals/validation.json.

PROMISING NEXT ACTIVE ATTACK
Recover hidden structure from H itself, rather than matching an exposed chain
basis. Read interaction_commutant_route.md: kernel A->[V,Q(A)] on Hermitian
onebodyoperators forquarticV; if maximalcommuting algebra withsimplejoint
spectrum, recoverdensitydensitybasis, thenrecognizegraph/localratio structure.
Needactualexactorcertifiedhiddenorbitalrotation experiment; NOTimplementedyet.
DimensionM aloneinsufficient (spinSU2nonabeliancommutant). NumPy eigvectors
arenotexactrotationcertificate. Couldstart rationalGivens-hidden ratiochain
M4/6/8, rootexisting orbitalrotationhelpers. PreservefullHamiltonianonlyinput.
Otherlivehypotheses: weakcouplinglinkedclusters withuniformprovedtail and
positive normalform; no suchfermionicbridgeproved. Notes primarysources and
conditionalcluster_expansion_certificate_bridge.md explain count/precision.

PROVENANCE / PITFALLS
* Remote executed archive results/lambda_runs/structural_scaling/source.tar.gz
  sha aa81856a565d3dbfd84b140cdda314734476cf9e9b25dc346bd3ca72ba0dfe4c.
  Sources remote /home/ubuntu/spectra-structural; venvspectra-venv/bin/python.
* Currentproofsources evolvedafterremotefreeze; final localproofhashes in
  structural_intervals/validation.json. Earlieragentdraftscontainerrors,
  superseded byrootcode; readreport andcurrent independentaudits.
* AlignedFCIuppergeneration active_space_reference_upper.py reconstructs
  h1 andERI fromfrozenrationalH: oppositespinquartic coefficient dividedby
  canonicalword SIGN, NO factor4. RebuildingSCF causedbasisphasesignmismatch;
  initialagentfailedformulaalsoforgot signandmultiplied4. Fixed and H4/H6
  energy expectation≈FCI checked1e-7. H8/H10top1000truncation contributeswidth.
  Useactive_space_ladder_references_aligned(_localized) only. FCIcostvalidation
  separate, no claimscalableupperstatediscovery. summary.json inlocalizedfolder
  lastinvocation listsH8/H10only; eachh4/h6receipt preserved, useperfolder.
* No production files changed. Old interruptedchainwork remainsseparatebelow.
* Condapython /opt/homebrew/Caskroom/miniconda/base/bin/python forSciPy/cvxpy,
  .venv-molecule/bin/python forPySCF; -S forstdlibproofs. Systempython3CAerror,
  condaSSLworks. NeverdisableTLS, exposeAPIkey, or tearhostsdownperbatch.

---

# Latest continuation — adaptive pricing completed

User said continue with next, following explicit authorization for many agents,
parallel Lambda experiments, API/CLI access, and two warm A10 hosts. No active
goal and no scheduled research automation. Do not contact archived tasks.

Read research/certificate_scaling/ADAPTIVE_RESULTS.md for this completed step.
32 remote experiments (24 original, 6 optimized, 2 ablations) exited0; all32
exported lower+upper intervals independently replayed using python -S. Five
local exploratory intervals also replayed. All130 downloaded file hashes and
four remote archive hashes matched. Six focused tests pass; production code
unchanged. Source snapshots/results: results/lambda_runs/adaptive_pricing/
and adaptive_pricing_v2/. Full numerical summary:
results/certificate_scaling/adaptive_intervals/campaign_summary.json.

Positive: dual-guided unequal-weight, width-four quadratic factors escape the
EXACT equal-weight-family ceiling on square H4. Best lower -3.471471331702718;
old family ceiling -4.459007306545172. New intervals: square .141082726552Ha,
rectangle .029386140487Ha, H6(best local cubic width8)3.486464894363Ha.
No new interval passes .0015Ha. Full quadratic SDP controls remain better
(.00503664, .000857081, .01302935Ha respectively), usually smaller as well.
No generic compactness, competitive advantage, FeMoco or superconductivity claim.
Full pair-cone controls optimal_inaccurate; valid primal bounds, no cone ceiling.

Measured pruning regression: square w4/full/no-prune81solves lower-3.47147 in72s;
half/no-prune81solves -3.47222 in69s; full/prune81solves -4.45901 in7.88s.
Half+prune201solves -4.38075 in23.4s. Pruning stays OFF by default. Row reduction
is valid for this Hermitian formulation with residual weights1/2, but dual
pricing trajectories can differ. Do not call the pruning speed a discovery win.
Both modes exactly replay complete polynomial residuals; full pricing maps and
matrices still built. No omitted-family optimality or bounded-cost theorem.

Next scientific step: test jointly adjustable small PSD blocks instead of fixed
rank-one directions, retaining historical constraints. Accuracy/size/discovery
cost must beat existing full-SDP controls across H4 geometries and H6. Sparse
block selection must ultimately avoid or account for complete pricing maps.
Do not resume tuning equal-weight pair count as if the proved obstruction vanished.

BOTH A10 HOSTS ACTIVE INTENTIONALLY, all current research processes completed:
129.159.32.200 (fab09ca8165b41d7bfc5649f413dda3f),
146.235.200.232 (b55288de9761445f8c589c21829f853c).
$1.29/h each, $2.58/h combined; keep warm, no per-batch teardown.
Use /opt/homebrew/Caskroom/miniconda/base/bin/python
research/certificate_scaling/lambda_pool.py status. System python3.13 has a
local CA-store error; conda Python verifies TLS successfully. Never disable TLS.
API key ~/.config/spectra/lambda-api-key is0600,outside repo. Never print it.
SSH userubuntu, ~/.ssh/id_ed25519. Remote env /home/ubuntu/spectra-venv/bin/python.
Executed sources /home/ubuntu/spectra-adaptive and spectra-adaptive-v2; archived
source hashes in ADAPTIVE_RESULTS.md. Jobs use CPU cores on A10 hosts, not GPU.
Prior direct discovery reference: research/certificate_scaling/DIRECT_RESULTS.md.
Older interrupted chain work below remains unfinished and separate.

---

# Latest completed campaign — certificate scaling, 2026-09-12

User explicitly authorized lots of subagents and Lambda experiments toward
compact efficiently discoverable fermionic certificates. Eight agents + root
completed a bounded campaign. Read research/certificate_scaling/RESULTS.md first.
No production code changed; earlier interruptedv20/v16 work below remains unfinished.

One Lambda A100 was launched at$1.99/hour under existing$2600 total authorization.
TEN jobs completed in42.85s batch.69 downloaded resultfile hashes matched.
37 returned lower certificates independently replayed via python -S locally.
Instance158.101.19.250 TERMINATED; provider UI No running instances verified.
Observed runningcost estimate~$0.29 (not invoice, unobserved interval excluded).
Source/results/lifecycle preserved results/lambda_runs/certificate_scaling/.
No ongoing GPU jobs. Never contact archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

REAL RESULT: round existing sparse factors todenominator100000, removezeros,
replayentirelower and exactintegerupper. H4square interval width.00127285521019Ha,
175912->119096bytes(32.30% smaller). Same frozenrule H4rectangle width
.000731085390761Ha,129032->112873bytes(12.52% smaller). Both<.0015Ha target.
Finite rationalSTO3G electronicHamiltonians only; no continuum/integralalgorithm
certification. Upperreplay70states, originaldiscoverycost not removed.
Six fullintervalcertificates inresults/certificate_scaling/intervals/summary.json.

NEGATIVE/MIXED: spectralranktruncation loseslargeaccuracy; idealLP repairhelps
fixedcompressedfactors butstillworseoriginal. Newrootfull_residual_replay.py exact
centeredonebodyresidualbound helpsunrepairedsome, negligibleafteridealrepair.
4modeadaptivepricing notbetterfixed3blockbaseline; all searchcostcounted.
60Lambda localitysettings(120disjoint/overlaprecords) numerical only;16 separate
integerclusterPSD exactcompositions locallyL4/16/32/64. Overlapstillnumerical.
SingleGPU H4endtoend5.57s vsCPU4.67s isnotcontrolledrepeatbenchmark andnoGPUwin.
Raw remoteheldout_gpu name MISLEADING: sameH4squaregeometry symmetrysplit,
coefficientrounded representation. Only rectangle isactualgeometrytransfer.

Earlyexplorationerrors caught and corrected: missingdaggerCARmap, wrongLP
solution slice, rationalfloatreportingerror, LDLpivotbug, invalidtwo-sidednumber
idealtheorem. Initialproxies/perblockreports notauthoritative; use RESULTS.md
and fresh exact receipts. Next actualtarget directfactor support/precision
search+constraintrepair with cost/errorcontrolled on growingactivespaces;
no claim generic efficientdiscovery or betterthanDePrince/Google demonstrated.

## Previous continuation follows

# Local continuation — 2026-09-12

## User steering and interrupted work — broader quantum-matter mission

USER STEERING supersedes the automatic habit of adding another pair of chain
constraints. User explicitly wants broad chemistry, FeMoco and room-temperature
superconductivity as the research north star. Ground energies alone are inadequate.
Pursue a reusable discovery/representation mechanism spanning energy/spectra,
forces/reaction free energies, and finite-temperature response/coherence. Do not
promise universal efficient solutions or material discovery. The pasted blueprint
contains unproved shortcuts: energy brackets cannot simply be differentiated into
force bounds; a single resolvent is not general causal finite-T response; replacing
energy with free energy does not automatically yield tractable certified entropy.
Relevant primary result checked: Fawzi/Fawzi/Scalet 2024, Nature Communications,
https://www.nature.com/articles/s41467-024-51592-3 . Certified thermal observable
hierarchies exist; fast convergence only established in restricted settings.
Do not resume incremental chain optimization as the whole strategy without tying
it to the new broader research mechanism. Preserve and finish verification of
already written code; do not claim interrupted work is a completed milestone.

INTERRUPTED v20/v16 work exists in results/marginal_graded_hubbard8/residual_joint/.
Five production/driver files changed, snapshots in source_before/original_sha256.json.
New experiments/marginal_residual_coherence.py implements last checkpoint's two
residual labels. ENERGYv20 adds field residual_coherence; FAMILYv16 imposes2 new
zero moments, cap211. Actual numerical rows202, energy coordinates201. New tests
89 focused PASSED457.35s and3 fractionfree PASSED. FULLSUITE NOT YET RUN.
No agents used. All observed numerical/energy jobs have finished; recheck live
processes before deciding any job terminal or restarting. No new jobs launched
since user requested strategic reassessment.

Equal500-evaluation-cap paired searches used same accepted v19 seed:
W0 control periodic-.6429067598293581 (rounded WORSE than seed), enlarged
-.6428975994922218; W1 control-.6605311930958622, enlarged-.6605207482896872.
All6 seed/control/enlarged matched energy replays ACCEPTED; seed/control v15 family
replays ACCEPTED. Numeric scope string in enlarged helper retained legacy199 label;
actual vectors and code have201 coordinates,77 offdiagonal. Preserve source hashes.
W0 enlarged initial40-round family exact proposal exists:200sources,
cap-.6428386735955448; production FAMILY replay/closure/fullRDM NOT run yet.
W1 exact family export FAILED inconsistent. Its ledger preserved; no acceptedv16cap.
W0/W1 family_resume scripts new202-row solver; do not relax existing4096-char gate.

held_out/{seed,control,enlarged} energy replays ACCEPTED. TargetU5,t1,V1/4,W-1/5.
Proposed periodic seed-.5223983238810757, control-.522556699829358,
enlarged-.5227293394922218: adverse transfer, must preserve signs.
Matched W0,W1 andheld_out without_residual proposals and congruence witnesses
created; exact ablation energy replays NOT launched. New comparison helpers
residual_joint_compare.py and residual_joint_transfer_compare.py written NOT run.
They require those ablation proofs. Newclosure/overlap helpers written NOT run.
No finalcollector, report or newprovenance for this unfinished milestone yet.
Prior607-file pure_coherence and116-file residual_coherence manifests verified
before production edits; historical snapshots required now. General goal ACTIVE.

## Latest accepted checkpoint — residual directions and a second frozen-recipe limit

This turn made VERIFIED PROGRESS; goal ACTIVE. All jobs TERMINAL. No production
changes, agents, GPU or paid resources. Do not contact archived task
01a08e58-e30c-7800-8df6-ab713cbd7787. Read research/marginal_residual_coherence_obstruction.md
and results/marginal_graded_hubbard8/residual_coherence/summary.json first. The prior
pure_coherence v19/v15 energies and full-family brackets below remain authoritative.

NEW EXACT DIRECTIONS (discovery only; no ENERGYv20):
- T0 from |358>+|409>: Y8 entries, T56 entries,8 changed bits. NormY<=1/4,
  normT and open boundary sum<=1/2. On selected W0 mixture, exact nonzero moment
  +.0001624428642127139.
- T1 from |103>+|358>: Y16 entries,T128 entries,2 changed bits: occupation-dependent
  same-spin hop across FOUR sites. NormY<=1/8,T/open<=1/4. W1 exact moment
  -.00013784965735813527. Do not use a generic >2-bit one-body exclusion argument.
Both pure offdiagonal operators reproduce independent accepted full-RDM mismatches;
removed diagonal moment zero from complete spin-word closure. Signed symmetries,
Hermiticity/spin, fermionic embeddings, cyclic cancellation and norm bounds exact.

INDEPENDENCE beyond all ENERGYv19/FAMILYv15:
residual_coherence_obstruction.py exact bounded row elimination finds
L0(M)=-16*M[413,1433], L1(M)=8*M[1382,1433]. Each annihilates all75 old offdiagonal
operators, including both pure_coherence terms. Sector(4,2) excludes actual fixed
half/charged projectors. Offdiagonal excludes all diagonal terms; neither selected
entry is a nearest-neighbor transition, excluding physical hopping profiles. New
minor [[0,1],[1,0]], determinant -1. Both directions independent of entire oldfamily.
Local positivity is intact; only stationary extension of these particular mixtures
is refuted. No generic representability conclusion.

EXACT FROZEN-RECIPE CAPS (all old v19 fields fixed; arbitrary REAL newcoefficients):
residual_coherence_fixed_limit.py <=20 pricing rounds, exact3-row fractionfree;
residual_coherence_fixed_limit_replay.py independently contracts physicalCAR,
ALL oldterms including pure_coherence, and fixed projector penalties. Trace1 and
both new moments exactlyzero imply lambda_min(A+g0T0+g1T1)<=Tr(rhoA). Fixed penalty
offset subtracted then divided by5. No fidelity inequality needed for fixedpenalties.
- W0 ceiling -.6429063846048632; gain above accepted seed at most
  1.992762124239207e-7/site. Three sources400/368/300 amplitudes,maxweight69chars.
- W1 ceiling -.6605360247761094; maximum gain7.146841714916643e-8/site.
  Three sources368/400/400 amplitudes,maxweight65chars.
These resolve only the frozen-recipe two-coefficient subproblem to the stated
bracket widths. Full v19 family gaps unchanged4.49049126240775e-5/.0001670212117076.
Exact proposals/replays under residual_coherence/<case>/fixed_limit/.

NUMERICS remain NONACCEPTING:
W0 proposednewcoeffs[7/1000000,1/8000],gain1.6e-7/site,9spectralevaluations.
W1[63/1000000,-3/500000],gain4e-8/site,14evaluations. 250budget, derivative and fresh
physical checks pass. No new energyPSD/schema or transferclaim. Prior adverse
whole-frozen-recipe heldout result remains unchanged.

VALIDATION/PROVENANCE:
40 focused tests PASS50.30s:22new +18old coherence/cap. Newtests include long-range
2-bit physical exclusion and nonzero old-pure-coherence contribution. Initial new
physical test fixture violated declared profile mean/reflection and was properly
refused; corrected to valid nonuniform profiles. Failed testsource/log/XML preserved
in residual_coherence/initial_fixture/. No proof source modified after receipts.
No fullsuite rerun: production unchanged, prior1128+102subtests remains historical.
Five newacceptedreceipts,355proofhashes,569construction/proofhashes audited. All607
prior pure_coherence provenancefiles unchanged. New116-file manifest VERIFIED;
central residual_coherence_obstruction matches summary; prior pure_coherence_energy
section matches unchanged prior summary. Report completed and opened.

NEXT: pursue a joint search with older fields free, or a broader coherence search,
to improve the full brackets; a frozen2-coordinate search is now quantitatively
limited to tiny gains. Do not infer useful joint gains solely from independence.
If adding productionv20, snapshot5 current v19 production/driver files beforeedits;
preserve both607 and116 manifests. Historical receipts would be291older+38v19+5new
=334. Do not overwrite old discovery sources/reports. Full-family attainment,
representability, broadertransfer and requested-accuracy scalability remain unproved.
Goal must stay ACTIVE, not complete or blocked.

## Latest accepted checkpoint — joint pure-coherence energy integration

This turn made VERIFIED PROGRESS; goal ACTIVE, not complete or blocked. All launched
jobs are TERMINAL. No agents, GPU or paid resources. Never contact archived task
01a08e58-e30c-7800-8df6-ab713cbd7787. Read research/marginal_pure_coherence_energy.md
and results/marginal_graded_hubbard8/pure_coherence/summary.json first.

ENERGYv19/FAMILYv15 now integrates both pure-coherence directions identified at the
previous checkpoint. New production module experiments/marginal_pure_coherence.py
reuses the old coherent-projector builder and removes diagonal entries. Labels:
346,409,1 and 314,614,1. Both enter actual PSD matrices; old versions refuse the new
field. FAMILYv15 imposes both exact zero moments and the full preceding hierarchy.
The conservative source cap is 209; independent numerical rows 200. Existing
4096-character weight gate, 94 blocks, 4096 states and local PSD cap200 unchanged.
Five modified old production/driver files preserved in pure_coherence/source_before.

SELECTED ACCEPTED EVIDENCE under pure_coherence/:
- W_zero/refined: periodic lower -.6429065838810756, million-site open lower
  -.6429090838810756, family ceiling -.6428616789684515, gap4.49049126240775e-5.
  183 positive sources, maximum rational weight1859 characters.
- W_plus_1/scaled: periodic lower -.6605360962445265, million-site open lower
  -.6605405962445265, family ceiling -.6603690750328189, gap.00016702121170760045.
  194 positive sources, maximum rational weight1864 characters.
Both selected energy certificates are byte-identical to the polished final/ CPs.
Lower improvements over v18: 5.53523367336e-6 and2.65871180638e-5/site. Physical
uppers unchanged -.6106763470511881/-.6184244823693281. Ceilings bound attainable
relaxation lower certificates, not physical ground-energy upper bounds.

COMPARISONS:
Previous v18 energy and v14 family proofs freshly replayed in previous_family/.
Signed separation from preceding full-family ceilings remains NEGATIVE:
-5.197799697763369e-5/-.000353425424980599. No whole-family separation proof.
Prior two-new-coefficient frozen-old-recipe caps reconstructed independently and
bound to matched physical context: joint lower exceeds them by
5.497541967328133e-6/2.641432111163255e-5. Older coefficients responding mattered.
Matched fixed-recipe removal of ONLY pure_coherence loses .00113616/.0010779/site.
Inherited directory name final/without_coherent actually means removing only the
new pure-coherence field, v19->v18. Comparison helper hardened to bind previous
accepted proof; earlier helper/four receipts preserved in comparison_before/*.before.

EXACT CLOSURE / REMAINING OBSTRUCTION:
Both selected mixtures close all120 spin-word and both new moments exactly, with
stationary classical order5 Markov extension. Independent full five-site density
matrices verify both original positive-projector differences are exactlyzero.
Residuals persist in full quantum overlap:
- W0: 29340 upper-triangle mismatches, zero diagonal; positive projector
  |358>+|409> detects +.0001624428642127139.
- W1: 30108 upper-triangle mismatches, zero diagonal; |103>+|358> detects
  -.00013784965735813527.
Only these particular symmetry-averaged mixtures are refuted, not all mixtures at
these energies or the ceilings. Local positivity intact. Fresh full test-space
census3960=120 diagonal+3840 offdiagonal. General quantum extension remains unproved.

FROZEN HELD-OUT TRANSFER:
Initial W0 recipe BEFORE polishing transferred to U5,t1,V1/4,W-1/5. All three
energy proofs accepted. New terms help +.00053568/site, but the complete frozen
new recipe is WORSE than the previous spin-word recipe by .00031310172234524/site.
Open lower new -.5223150350391923; ablated -.5228507150391922; previous
-.5220019333168471. Physical upper unchanged -.4885616802989547. Preserve both
signs. No generic robustness or molecular/geometry/filling transfer claim.

NUMERICAL ATTEMPTS / REFUSALS:
Joint199 coefficients:120 diagonal,75 offdiagonal,2 hopping profiles,2 penalties.
Initial searches460/500 full-spectrum evaluations, polish461/500; no convergence
claim. Initial40 two-eigenvector family pricing yielded accepted ceilings
-.642830074068934/-.6602582074343529, retained in final/. Pruned ledgers retain
all accepted physical sources plus4096 anchors. Refined80 one-eigenvector rounds:
W0 exact export accepted412.38s; W1 exact basis inconsistent, diagnostic370.11s
preserved, no accepted ceiling from that attempt. NEW discovery clone
pure_coherence_family_scaled.py uses ZERO extra pricing, 1000 row scaling and
highs-ds/no presolve on the refined5064-candidate W1 ledger. Exact export accepted
190.52s; production physical replay, closure and fullRDM passed unchanged gates.
Initial scaled energy CLI call missed required arguments; failed log preserved,
corrected invocation succeeded. No frozen production/helper source overwritten.

VALIDATION / PROVENANCE:
74 focused integration tests,3 fraction-free tests pass. Full regression PASSED
1128 tests and102 subtests in1442.11s; existing calibration return-value warning.
No production changes after full-suite collection. Copied dimension test initially
retained old row boundary; corrected to200 accepted/201 refused, failure log kept.
Collector audits38 current accepted receipts /1591 proof hashes,1856 complete
construction/proof source entries. All607 new provenance files VERIFIED. Historical
291 receipts:9792 unchanged hashes+455 preserved-source hashes. Previous spin_word
manifest543 unchanged+5 snapshots (548); spin_coherence111+5 (116). Central new
pure_coherence_energy section matches summary; preceding two sections match their
unchanged summaries. Report and provenance completed; all jobs terminal.

NEXT substantive work: address remaining independent offdiagonal obstructions or
improve the family numerical bracket while preserving exact gates. First assess
new residual projectors and frozen-recipe gain limits, as at the prior coherence
stage; do not infer that merely adding a residual guarantees useful joint gain.
Full-family attainment, general representability, robust transfer and requested-
accuracy scalability remain unproved. Goal must stay ACTIVE.

## Latest accepted checkpoint — compact coherence directions and frozen-recipe caps

This turn made VERIFIED PROGRESS; goal ACTIVE. All launched jobs TERMINAL. No
production changes, agents, GPU or paidresources. Never contact archived task
01a08e58-e30c-7800-8df6-ab713cbd7787. Read research/marginal_spin_coherence_obstruction.md
and results/marginal_graded_hubbard8/spin_coherence/summary.json first. Previous
spin_word checkpoint below remains authoritative for all energy/full-family bounds.

NEW EXACT DIAGNOSTICS (discovery only; no ENERGYv19 yet):
spin_coherence_obstruction.py reuses full_overlap_telescope.build and removes diagonal
from the two newest positiveprojector telescopes. Complete120diagonal closure makes
removed moment exactlyzero in selected spin_word mixtures; purecoherence still
matches independently reconstructed fullRDM mismatch.
T0 source |346>+|409>:5siteY8nonzeros,6siteT60nonzeros,4changedfermionbits,
exactexpectation W0-.00021964912486317672; occupation-dependentspinexchange.
T1 source |314>+|614>:Y16/T120nonzeros,6changedbits,
exactexpectation W1-.00027092120898647947.
BothY row-sum normbound1/4; bothT and open translatedboundarysum normbound1/2.
ExactHermiticity/spin/symmetries/fermionictranslation/6sitecyclicsum verified.
Originalpositiveprojectorsdetectstationarity, not negative localpositivity.

INDEPENDENCE beyondfull ENERGYv18/FAMILYv14:
L0(M)=M[1370,1433]-M[350,413], L1(M)=M[1337,1637].
Both annihilate all73 oldnondiagonaloperators. Entriessector(4,2) excludesactual
fixedhalf/chargedprojectors;4/6bits excludesphysicalonebodyhopping;offdiagonal
excludesallold/newdiagonals. New2x2minor diag(1/8,1/16),det1/128.
ThusT0/T1 are independentmodulofulloldfamily. T0 hasno new supportentry alone;
its functional cancels two equal oldspin0,3 entries. Don't reuseold8-bit-only
supportargument. Crossmoments W0[-.000219649124863,-.0000641865373040],
W1[+.0000670478328833,-.0002709212089865].
Exactreceipts spin_coherence/<case>/telescope_replay.json andsupport_replay.json.

NEW EXACT FIXED-RECIPE LIMITATION:
spin_coherence_fixed_numeric.py two-coefficientsearch keepsALLoldfields frozen.
W0 proposed[0,0] zero gain;W1[-93/500000,179/500000], proposed1.4e-7/sitegain.
These are NONACCEPTING proposals, no new energyPSD orprotocol acceptance.
spin_coherence_fixed_limit.py <=20round three-rowpricing exactpositivebasis export;
spin_coherence_fixed_limit_replay.py independentstdlib physical CAR expectation
thenexactcap for ARBITRARYREALnewcoefs, with everyoldfield includingpenaltiesfixed.
Mixturetrace1 andbothnewmomzero give lambda_min(A+g0T0+g1T1)<=Tr(rhoA).
Subtractfixedpenaltyoffset/divide5. NO fidelityinequality required becausepenalties
arefixed. These caps DO NOT cover jointreoptimization ofoldfields orfullnewfamily.
- W0cap-.6429120814230429, maxgainoveracceptedseed3.7691706031867335e-8/site;
  ONE positive368-amplitude source in(3,3), exactweight1.
- W1cap-.66056251056563819, maxgain1.7279695216745168e-7/site;
  THREEpositive sources400/400/300amplitudes, sectors(3,3)/(3,3)/(3,4),
  maxweight54chars. Bothnewmoments exactlyzero. Accepted fixed_limit/fixed_family_replay.json.
Thereforevaryingonlynew2coefficients cannotyield largeimprovements; meaningfulgain
requires oldercoefficientsreopt and/orfurtherdirections. Fulloldfamilygapsunchanged
5.7513230651e-5/.000380012543044; fullquantumextension stillunproved.

VALIDATION/PROVENANCE:
30focusedoverlap/coherence tests +9fixed-limit tests PASSED. Includes18newtests,
independentpartialtrace/norm, diagonalremoval, directvacuumenergy, amplitude scaling,
nonzero moment andwrongscope/malformedmixture refusals. No fullsuite rerun because
no productionchanges; preceding1074tests+102subtests remainshistoricalfullpass.
Initial newfixedreplay misread _actions dictaslist, TypeError and2focusedfailures;
correctedexplicit4096stateindexing. Initialsource/logs preserved source_attempts/;
corrected9tests andbothactualphysicalreplays accepted. No failedreceipt promoted.
spin_coherence_collect.py audits5newacceptedreceipts/350proofsourceentries,
560allconstruction/proofentries; all548prior spin_wordprovenancefilesunchanged.
New116fileprovenance VERIFIED; centralspin_coherence_obstruction matchesnewsummary,
centralspin_word_energy matchesunchangedprior summary. Alljobsfinished.
No newenergybound or heldouttransfer result thisturn.

NEXT substantive step: integrate both purecoherence directions intojointenergy and
family optimization (wouldadd2to197energycoordinates and2to198independentfamilyrows;
existingcap207 growsconservatively2). Preserveoldversionrefusals, exactweightgate,
fullspinwordclosure andPSDblockcaps. Snapshot existing5v18production/driverfiles
beforeediting to preserve548/116manifests andhistorical291receipts (254historical
beforev18+32spin_word+5newdiagnostics). Reusepurebuildfromoldcoherent builder rather
than duplicatingalgebra. NewfieldmustbepartofactualPSDreplay; fixedrecipeproposals
areunsupportedandcannotmerelybe relabelled. Oldercoefficients must be free torespond;
a2coef frozen search is now quantitativelycapped and not the next useful step.
Otheroffdiagonalconstraints maymatter; genericrepresentability, jointnumerical
attainment, broadertransfer andrequestedaccuracyscalability remainunproved.
Goal mustremainACTIVE, notcompleteorblocked.

## Latest accepted checkpoint — full spin-word closure and stronger energies

This turn made VERIFIED PROGRESS. Goal ACTIVE, not complete/blocked. All launched
jobs TERMINAL. No agents, GPU or paid resources. Full regression PASSED1074tests
and102subtests in1352.97s; same existing calibration warning.63focused+3fractionfree
passed and included. No production changes after full collection. Never contact
archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

Read research/marginal_spin_word_energy.md and
results/marginal_graded_hubbard8/spin_word/summary.json first. New ENERGYv18/FAMILYv14
implements COMPLETE120 diagonal five-site spin-word directions: PH/spin even,
reflection odd. Preserves full earlier hierarchy, refusals,94blocks/4096states,
maxlocalPSD200 and4096-character weight bound. Newconservative sourcecap207; actual
independent numerical198rows (old9 selectedshapes redundant in120basis).
Fullspinword module experiments/marginal_spin_word_telescope.py. Five changed old
sources preserved in spin_word/source_before. Historical coherent sources/reports
not overwritten. General repr/scalability remains unproved.

SELECTED EVIDENCE under results/marginal_graded_hubbard8/spin_word/:
W_zero/scaled/: periodic lower-.64291211911474899, million-site open-.64291461911474901,
family ceiling-.64285460588409793, gap5.751323065099369e-5,192positive sources,
maxweight2103chars. Lower improvement overv17 1.353189346852e-5/site.
W_plus_1/final/: periodic-.66056268336259027, open-.66056718336259035,
family ceiling-.66018267081954596, gap.00038001254304439894,198positive sources,
maxweight2768chars. Lower improvement4.357950876808e-5/site.
Physical uppers unchanged-.6106763470511881/-.6184244823693281. Familyceilings bound
attainable relaxation LOWER certificates, not physicalgroundenergy uppers.
Exact previous-family comparisons nonseparating (-2.1600537065e-5/-1.65402962809e-4).
Matched fixedrecipe removing ONLYspinword field loses.03243578/.01586254/site.
These large fixedrecipe gains do not establish reoptimizedwholefamily separation.
Precedingenergy+familyreplays fresh in previous_family/. Matchedablation proofs in
final/ bound to selectedscaledenergy by exactcertificate byteidentity.

Both selected mixtures have exactlymatching full1024 five-site spin-word diagonal
laws,120zero moments and classical stationary order5 Markov extensions on
empty/up/down/double. W0positive Markovstates960/flows2566;W1states1024/flows3228.
Proof: difference of prefix/suffixdiagonals lies in120-dimensional symmetryspace;
all120telescopemoments zero forceit zero. Independent4096-integerlawconstruction
confirms. No offdiagonal coherences or fixedglobalparticle number extended.
FullRDMdiagnostics accepted, zero diagonal residuals but QUANTUMextensionrefuted:
W0 30136uppermismatches; positiveprojector v=|346>+|409>, mismatch-.00021964912486317672.
W1 30324uppermismatches; v=|314>+|614>, mismatch-.00027092120898647947.
Local positivity intact; refutes onlythese particularsymmetryaveragedmixtures.
Exactfullrealoverlap census3960=120diag+3840offdiag, not scalability proof.

Frozen INITIAL W0 recipe beforepolish at U5,t1,V1/4,W-1/5 in held_out/:
withspinword open-.52200193331684708; without-.53905959331684705;
previouscoherent frozen-.52245353505179537. Signed newfield contribution+.01705766;
wholefrozenrecipe improvement+.00045160173494836. All3energyproofs+exactfrozenfield
comparison accepted. Prior coherent2term harmfultransfer remainshistoricalevidence.
Inherited directories without_coherent actually remove ONLYspin_word_telescope;
comparisons verifyfieldidentity. Physicalupper-.4885616802989547unchanged.

NUMERICAL LIMITATIONS PRESERVED:
197energycoordinates,500evalcap. InitialW0/W1 410/500eval;polish417/500.
InitialW1 family40twoeig accepted198sources. InitialW0 LPfailed HiGHSUnknown;
separate spin_word_family_trust.py boundedfallback recovered accepted196sources
in retry/ (copied tofinal/). All originalfailurelogs retained. InitialW0cap
-.6428006209165888, gap.00011149819816011205.
Both refined80oneeig ledgers preserved acceptedphysicalsources+determinantanchors.
W0/W1 numericalcaps-.6428546060245608/-.6603172578263036 BUT exactexportsREFUSED
inconsistent189/195-column bases. Tinyfloatingresiduals never accepted. Those
refined dirs have acceptedenergyonly, no familyproof/proposal. Newzero-pricing
spin_word_family_scaled.py uses allconstraints scaled1000, highs-ds nopresolve,
originalexact198-row solver andunchangedweightgate. W0exact192-source reconstruction
and physicalfamily+closure+fullRDM passed (selected scaled/). W1 exact196-column
reconstruction stillINCONSISTENT, nofamilyproposal; retaininitialacceptedfinal/.
No enlargedfamilyconvergence/attainment claim. Further shrinkage requires addressing
exactbasisconditioning, not assuming numericalLPoptimalstatus is acceptance.
spin_word_lp_scale_probe.json records nonaccepting scale1/1000/1e6diagnostics.

spin_word_collect.py accepted32currentreceipts/1292sourcehashentries;254historical
receipts with8226unchanged/379snapshot entries. Prior481filemanifest accountedfor
476unchanged/5snapshot. New548fileprovenance VERIFIED; centralspin_word_energy
matches summary and precedingcoherentsection/summary unchanged. Newreport and
central validation section written. No currentfailedproposal promoted toproof.

Next substantive work: offdiagonalcoherenceconsistency beyondexisting73directions,
compact newprojectortelescopes fromresidualwitnesses, and/or stabilize exactfamily
reconstruction (remaining gaps above). Full3840offdiagonal constraintcost isunproved.
All diagonal occupationconsistency isnow closed forthese symmetryaveragedfamilies;
no reason tokeep addingdiagonalwords. Generalquantumrepresentability, numerical
attainment, broaderphysicaltransfer andrequested-accuracyscalability remainunproved.
Do not markgoalcomplete/blocked.

## Latest accepted checkpoint — coherent-projector energy and transfer limitation

This goal turn made VERIFIED PROGRESS. Goal ACTIVE. Production ENERGYv17/FAMILYv13
and current exact proofs are implemented. FULL REGRESSION PASSED:1039tests and
102subtests in1191.30s, same existingcalibrationreturnvaluewarning. All launched
jobs are TERMINAL. No agents, GPU or paid resources. Prior archived task unchanged.

New results root: results/marginal_graded_hubbard8/coherent_projector/.
Selected proofs in <case>/final/. New production module
experiments/marginal_coherent_projector_telescope.py fixes two labels
358,601,-1 and346,613,1. Five existing sources snapshotted before edits in
coherent_projector/source_before/. Old versions/refusals/sourcecaps preserved;
new139-source cap and two exact zero moments require the full prior hierarchy.
Local94blocks/4096states/maxPSD200 unchanged; weightlimit4096characters unchanged.

W_zero: periodic -0.64292565100821752; open -0.64292815100821743;
family cap-0.64289051857768398, gap3.5132430533492614e-05;139positive sources.
Lower improvement6.1237590234e-06; strict previous-family separation FALSE.
Matched fixed-recipe removal of newterms loses0.0005827/site.
Both newmoments exactlyzero but fullRDM stillfails; spin-word diagonalTV0.00083285600913601614.
W_plus_1: periodic -0.66060626287135837; open -0.66061076287135845;
family cap-0.66039728039978174, gap0.00020898247157668336;139positive sources.
Lower improvement2.57706019341e-05; strict previous-family separation FALSE.
Matched fixed-recipe removal of newterms loses0.00020318/site.
Both newmoments exactlyzero but fullRDM stillfails; spin-word diagonalTV0.0019415404864194311.

W0 residual full-overlap v=|102>+|153>, mismatch-.00024514857339673666,30760upper differences/616diagonals.
W1 residual v=|103>+|358>, mismatch-.0002132191309802359,31000upper/664diagonals.
Coarse charge-word laws allagree exactly; full spin-resolved wordlaws fail.
Diagonal positiveprojectors308/332states detect the TV above; norm<=1 telescopes.
Exact signed-orbit census: real spin-conserving5site Hermitian dimension32264;
PH-even/spinflip-even/reflection-odd3960=120diagonal+3840offdiagonal directions.
This dimension count is not a scalability proof. Full density/diagonal probes
and independentcoherentclosures accepted using stdlib exactinteger/rationals.

Frozen transfer uses INITIAL W0 beforepolish: held_out/. U5,t1,V1/4,W-1/5.
Withnewterms open-.5224535350517954; without-.5216497950517954;
previousfrozen three-spectator-.5225841142675748. Signed newterm contribution
NEGATIVE -.00080374: newterms HARM this frozenheldout recipe. Entire newrecipe
improvesprevious by.0001305792157794 due jointlyadaptedothercoefficients.
All3energyreplays and frozencomparison accepted. Initial matched-targetdriver
correctlyrefused missingheldout filterrecipe; logsretained. Correct dedicated
coherent_projector_target_replay.py transferred fixedW0filter and evaluated
actualheldoutHamiltonian; didnotweaken matched-sourcegate.

Numerics:138coefficients,500evalcap. InitialW0/W1 412/500;polish411/500.
Initialfamily40two-eigenvectorrounds;refined80one-eigenvectorrounds. Finalledgers
11838/12789vectors;139sources,maxweight2756/2906chars. Fraction-free139rowhelper
bounded4096-bit clearedinput/50000-bit pivots and exactoriginalequations.
W0/W1 refineddurations485.61/468.72s,notcontrolledbenchmarks. Last reduced
eigenvalues remainnegative; no convergenceclaim. Bothpreviousfamilycapsfreshly
replayed undercurrentcode and signedcomparisons accepted, nonseparating.

81focusedtests and3fractionfreetests passed; both included in full1039suite.
coherent_projector_collect.py accepted29currentreceipts/1119sourcehashentries;
225historicalreceipts with7157unchanged/329preservedsnapshot entries. Newprovenance
481files. Read research/marginal_coherent_projector_energy.md and coherent_projector/
summary.json first. Central coherent_projector_energy section added; oldreports
and summaries preserved. No productionchanges afterfullcollection. Final481file
hashes verified, currentcentralmatches and priorsummariesunchanged; alljobs terminal.

Next substantive work: broaden stationary spin-word constraints (120 diagonal
directions) and/or full3960-dimensional real overlap space rather than assuming
new2 scalarconditions finishquantumconsistency. Tightenfamilygaps; transfer
robustness needswork because signedablation isadverse. Generalrepresentability,
genericmolecular/longrange/higherDtransfer and requestedaccuracyscalability remain
unproved. Do not markgoalcomplete orblocked.

## Latest accepted checkpoint — full quantum overlap obstruction

Read `research/marginal_full_overlap_obstruction.md` and
`results/marginal_graded_hubbard8/full_overlap/summary.json` first.
This turn made VERIFIED PROGRESS. Goal ACTIVE, not complete/blocked.
All launched jobs terminal. No production edits, agents, GPU or paid resources.
Never contact archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

The following preceding three-spectator checkpoint remains the latest ENERGYv16/
FAMILYv12 energy evidence. Its prior statement that full-RDM diagnostics were
not implemented is now superseded. Both selected137-source local mixtures FAIL
full five-site overlap after eight-image PH/spin/reflection symmetrization.
Exact positive projectors v v*/||v||² detect:
- W0 v=|358>-|601>, norm²2, sector(3,2), left.03455734832389095,
  right.034272965597014836, mismatch+.0002843827268761191.30952nonzero
  upper difference entries,648diagonals,4709-bit common denominator.
- W1 v=|346>+|613>, norm²2, sector(3,2), left.00782211122115541,
  right.00807266225122646, mismatch-.0002505510300710504.31000upper,
  664diagonals,4926-bit denominator.
All probabilities valid: overlap inconsistency, NOT negative local positivity.
Spin-resolved diagonal mismatch does not contradict old charge-Markov closure.
This refutes a stationary quantum extension of THESE particular symmetry-
averaged mixtures; no general no-extension or numerical-optimum conclusion.

Discovery full_overlap_density.py reconstructs complete contiguous five-site
partial traces with bounded integer arithmetic, verifies source orbit closure,
unit traces and spin-sector support, records hashes and exact sparse witnesses.
full_overlap_telescope.py makes A=PH/spin average(P),Y=(A-RAR*)/2,T=Yleft-Yright.
EachY16nonzeros,denominator16;T W0 116nonzeros/W1 124nonzeros. ExactHermitian,
spin-preserving,PH/spin/reflection-invariant, fermionic left-to-right translation
and six-site periodic sumzero. Analytic||Y||<=1/2,||T||<=1,openboundary<=1.
Five-sitePH²=-I is projective, tested explicitly. Direct exact original137source
expectation equals independently reduced symmetry-averaged mismatch.
Both telescope replays and support replay run python-S (stdlib only).

full_overlap_support_replay.py checks all71prior nondiagonal telescope actions
on4096states and actual fixed half/charged projector sectors. NewT0entry
(1382,1625)=-1/8;newT1entry(1370,1637)=+1/8, bothsector(4,2),8changedbits.
Twooperator minor determinant-1/64 proves TWO INDEPENDENT NEW DIRECTIONS modulo
fulloldaffinefamily. Four-spin coherences absent from old one-hop and spin/pair
exchange directions. Independence does not guarantee energy improvement.

21focusedtests passed. Prior982full+102subtests and3laterfraction-free tests
remain historical; no fullsuite rerun thisturn and no production change.
Five new acceptedreceipts in full_overlap: twofull_overlap_replay,twotelescope,
support_replay. Intermediate projector_witness.json staysacceptedfalse proposal.
full_overlap_collect.py audits allfive/sourcehashes and prior479-fileprovenance;
newsummary/provenance and centralfull_overlap_obstruction section added.
Allpriorreportsandmanifests preserved. No newenergy/familycap/transferclaim.

NEXT: integrate these TWO exact sparse projector telescopes in production energy
and family paths (newversion needed, oldrefusals preserved); searchcoefficients,
fresh exactPSDandmomentreplay, assessimprovement. Productionchangeswillstale
receipt hashes: snapshot oldsources first. Remaining currentfamilygaps W0
2.7128835223944408e-5/W1 .0005974181820038868 unchanged. Fullrepresentability,
numericalattainment, genericmolecular/longrange/higherDtransfer and requested-
accuracyscalability remainunproved. Donotmarkgoalcompleteorblocked.

## Latest accepted checkpoint — three-spectator energy, caps and transfer

Read `research/marginal_three_spectator_energy.md` and
`results/marginal_graded_hubbard8/three_spectator/summary.json` first.
This turn made VERIFIED PROGRESS. Goal ACTIVE, not complete/blocked.
All launched numerical, replay and test jobs are terminal. No agents, GPU or
paid resources used. Never contact archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

LATEST BOTH CASES: `three_spectator/<case>/final/` under the graded result root.
ENERGYv16/FAMILYv12 integrate all18 prior three-spectator obstructions.
- W0 exact periodic lower -.64293177476724088; open -.64293427476724088;
  physical upper -.6106763470511881. New family cap -.6429046459320169,
  gap2.7128835223944408e-5;137positive physical sources. Lower improvement
  over preceding certificate1.2123576518e-5.
- W1 exact periodic lower -.66063203347329248; open -.66063653347329248;
  physical upper -.6184244823693281. New family cap -.6600346152912886,
  gap.0005974181820038868;137positive physical sources. Lower improvement
  7.916114251224e-5. Neither enlarged numerical limit exactly resolved.
Family ceilings are NOT physical energy uppers. All final energy/family and
independent18three/30two/14one/4pair/4spin/coherent closures accepted. All these
moments zero; averaged charge laws have classical Markov extensions. Full
quantum extension unproved. No new violations tested in these newest mixtures.

W0 STRICT preceding-family separation is proved: sharper old two-spectator
cap-.6429325359567474 (111sources) freshly replayed in W_zero/previous_refined/.
New lower exceeds it by7.611895064903942e-7. Matching target, size, sources,
overlap ceilings, ratio and sparse span checked. W1 whole-family separation
NOT proved: signed excess vs old accepted cap-.6605062601234778 is
-.00012577334981470596. W1 previous_primal/ has only a nonaccepting inaccurate
conic objective~- .66053622; it was not promoted to a new cap. Do not confuse
improvement over an old accepted lower with improvement over its whole family.

Frozen transfer uses INITIAL accepted W0root, not polished/final recipe.
U5,t1,V1/4,W-1/5 on same half-filled chain: open lower-.5225841142675748;
upper-.4885616802989547. Removing only18newterms gives-.5254907142675748:
exact fixed-recipe loss14533/5000000=.0029066. Previous frozen two-spectator
recipe lower-.5226432283437589; improvement5.911407618408e-5. All3energy
replays and frozen comparison accepted under held_out/. No new geometry,
longrange or generic molecular transfer. Fixed-recipe ablation is not a cap
for a reoptimized preceding family.

Production: experiments/marginal_three_spectator_hopping.py has18canonical
labels i,j,k,l,m,p,q,r and exact combined actions. ENERGYv16 adds projected
matrix before PSD, old versions refuse new field. FAMILYv12 enforces18exact
zero moments and full preceding hierarchy. Fixed coefficient fields forbidden
in family inputs. Old modes/caps/refusals preserved; new sourcecap137 versus
119 in preceding mode. Rational weight limit4096characters unchanged. Local
94blocks/4096Fockstates/maxPSD200 and openingcost unchanged. Changed5existing
sources preserved in three_spectator/source_before/original_sha256.json.

Discovery: three_spectator_thermal_numeric.py uses136coefficients and500spectral
evaluation cap. Both initialruns500eval. W0polish390eval, W1polish500eval.
Selected polished lower proofs freshly generated congruence witnesses and
passed exact replay. Initial stale-witness attempt was interrupted, never
accepted; obsolete witness and interrupted log retained.

First40round candidate-feasible box runs exported nothing due missinghashlib
import; failedsource copies/logs retained, importfixed. Two global runs sharing
the bug were interrupted before export and restarted. Correct global40round
two-eigenvector runs produced W0 10514candidates, cap-.6428468597850456,
136sources,4091maxweightchars, accepted in W0global_final/. W1 10794candidates,
numerical cap-.6591750977671775,4217maxweightchars: PRODUCTION REFUSED weight
limit, recorded in W1global/weight_limit_refusal.json. Never accepted as a cap.
Global run durations387.24/416.54s, with costly rational reconstruction.

three_spectator_family_resume.py retains global ledgers and adds40candidate-
feasible trust-pricing rounds, two eigenvectors/block. New vectors use integer
scale10^6, old vectors retained. Outputs10869/11331candidates, selected137-source
caps and maxweight2834/2964characters. Run225.77/206.64s, not controlled speedup
benchmarks. New three_spectator_fraction_free.py solves at most137rows with
exact denominator clearing/divisions and original-equation check;4096-bit
cleared inputs/50000-bit intermediate bound. Existing85-row helper unchanged.
Last reduced eigenvalues stillnegative; no numerical convergence assertion.
Newest ledgers under13000 input bound, but do not blindly repeat pricing.

W0 preceding-cap refinement: existing two_spectator_prune_ledger.py pruned
9276old candidates to4215 preserving119sources+4096anchors;40trust rounds
produced4537candidates and accepted111-source cap above. W1 old-family conic
attempt is nonaccepting and did not establish required separation.

Validation: full982tests+102subtests passed in1099.08s, same calibration
return-value warning. Focused58passed in185.57s. Three later fraction-free
discovery tests passed separately in.40s; no full985-test claim. No production
change after full suite collection. Current38proofreceipts/1368sourcehashentries
audited;182historical receipts with5684unchanged and265snapshot hashentries.
Four snapshot directories resolve historical changes: pair_transfer/source_before,
pair_transfer/family_source_before, two_spectator/source_before, three_spectator/
source_before. New provenance479files, constructionhashentries1478. Historical
receipts are not fresh current-code replays. All prior reports/summaries kept.
Central marginal_final_validation.json has three_spectator_energy section.

Next substantive questions: tighten enlarged-family gaps, especially W1; test
full five-site reduced-density-matrix overlap agreement to locate quantum
consistency conditions beyond the scalar moments now closed. Such a full-RDM
diagnostic is NOT IMPLEMENTED. No claim of complete overlap/positivity closure,
general representability, generic molecular/longrange/higherD transfer or cost
at requested accuracy. Continue locally; goal remains active.

## Latest accepted checkpoint — W1 gap refinement and eighteen new obstructions

Read `research/marginal_two_spectator_refinement.md` and
`results/marginal_graded_hubbard8/two_spectator/refinement_summary.json` first.
This turn made VERIFIED PROGRESS. Goal ACTIVE, not complete/blocked.
All launched numerical/replay/test handles are terminal. No agents, GPU or paid
resources used. Never contact archived task01a08e58-e30c-7800-8df6-ab713cbd7787.

LATEST W1: `two_spectator/W_plus_1/enriched_final/`. W0 remains W_zero/final/.
W1 unchanged periodic lower -.6607111946158047; new family ceiling
-.6605062601234778; exact gap .00020493449232694596 (96.30624% reduction
from .005548129720482485).118 positive physical sources; max weight2534chars.
Open lower -.6607156946158047 and physical upper -.6184244823693281 unchanged.
W0 gap remains1.9009152447823635e-5. Neither limit exactly attained; W1 remains
unresolved. Family ceiling is not a physical energy upper. No new energy lower
or coupling-transfer run this turn; prior exact frozen transfer preserved.

New candidate-feasible trust pricing replaces coefficient clipping with an LP
that enforces all current candidate inequalities inside radius.05 dual boxes,
including both fidelity signs. Final production proof uses no trust boxes.
two_spectator_prune_ledger.py preserved every accepted118 source and4096anchors,
pruning14630 candidates to4211 (overlap between those sets). trust_seed/ immutable.
two_spectator_family_trust.py: first40 rounds2eigenvectors/block ->4749candidates,
accepted117-source cap-.658290338631094 in trust/. Next80rounds1eigenvector/block
->5445candidates, cap-.6602065679378636, consolidated/replayed in trust_final/.
Last reduced eigenvalues -.04939/-.01891; no pricing convergence. Both runs
use all94blocks. CPU durations64.45/120.29s, not controlled benchmarks.

Coherent subspace work: two_spectator_primal_subspace.py builds a primal SDP
over selected vector spans plusanchors,7blocks maxdimension39, one45solver-sec
run. optimal_inaccurate objective-.6604537282. two_spectator_primal_enrich.py
does4conic rounds with full94block pricing, capped32activeblocks/48dimension,
30solver-sec/80iterations per round; actualmaxdimension45. All statuses
optimal_inaccurate, finalobjective-.6605066556, PSD min~-2.25e-8/equality
residual~4.53e-8; some enrichments refused by32block cap. Global reduced
eigenvalues large negative, no global conic optimality claim.
two_spectator_family_atom_resume.py rounds up to16 positive conic atoms with
25variants each, keeps oldledger, adds<=400physicalvectors. Enriched run added
400 ->5845candidates, exact119-column reconstruction (118physical+1slack)
accepted. Its proposal in enriched_atoms/ supplied enriched_final/. The latter
uses unchanged polished lower and matched range-two profile. All10final
energy/family/separation/closure/probe receipts accepted.

Exact numerical-conversion issue resolved separately: primal_atoms/ kept5720
candidates but LP support118 could not span119exactrows; rectangular solver
correctly refused. one_column_physical_completion.py computes exactleftnull,
nominates physical columns (or valid fidelity slacks), tries<=4, requires exact
119equations and nonnegative weights. QR only nominates independent rows.
two_spectator_family_complete.py succeeded on firstanchor(index13,mass~2.833e-12),
yielding118physicalsources and cap-.6604538639170363. primal_completed/ energy
andfamily receipts freshly accepted. This is weaker than selected enriched cap.
Originalfailedattempt kept. Diagnostic basis_size118 is BEFORE completion;
physical_completion records addedcolumn; finalbasis119. No tolerance acceptance.

NEW obstruction: three_spectator_overlap_probe.py (18labels i,j,k,l,m,p,q,r).
All three remaining five-site positions are spectators.40PH-allowed products,
4self-reflected vanish,36remaining form18 antisymmetric pairs. T=Yleft-Yright.
FullCAR4096states/Hermiticity/PH/reflection/spinflip and periodic cancellation
checked; independentbit-swap tests cover every action. Norm bound8. Exact
restricted two-bit hopping Gram (N<=4) oldrank47 ->65: all18independent modulo
currentfamily, since fixedprojectors onlyN5/6/7 and otheroperators noentries.
All18momentsNONZERO in W0final and W1enriched_final mixtures. Largestmoments
W0 .00011186449945432143, W1 -.0005222541791278473; norm-normalized violations
1.3983062431790179e-5 /6.528177239098091e-5. Existing30two/14one/4pair/4spin/
coherent moments stillzero. ClassicalchargeMarkov extension stillvalid; no
stationary quantum extension of these particular averaged mixtures exists.
This does not invalidate their family caps. Three-spectator ENERGY INTEGRATION
NOT DONE. Next substantive direction: integrate/exploit these18conditions and
continue to addressW1remainingfamilygap; not blind unbounded pricing repetition.

Production unchanged ENERGYv15/FAMILYv11/sourcecap119. Previous945tests+102
subtests stillapply;7newfocuseddiscovery tests passed separately in2.58s.
Tests cover independent three-spectator actions/classes plus exactcompletion
including below-float residuals and negative/dependent refusals. No newfullsuite
wasrun. New26proofreceipts/852sourcehashentries audited, preceding26current
receipts/876entries re-audited; prior381-file provenance unchanged. New
refinement_provenance.json has198files, constructionhashentries856. Central
marginal_final_validation.json has two_spectator_refinement section.
Oldsummary/provenance/reports are preserved; newcollector is
two_spectator_refinement_collect.py. Allnew discovery sources under discovery/.

General quantum representability, generic molecular/longrange/higherD transfer
and requested-accuracy scalability remain unproved. Continue locally.

## Latest accepted checkpoint — two-spectator energy, family caps and transfer

Read `research/marginal_two_spectator_energy.md` and
`results/marginal_graded_hubbard8/two_spectator/summary.json` first.
This goal turn made VERIFIED PROGRESS. Goal ACTIVE, not complete/blocked.
All numerical, exact replay and test handles are terminal; process scan found
no remaining research Python jobs. No agents, GPU or paid resources used.
Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

Latest authoritative directories: `two_spectator/<case>/final/` under
`results/marginal_graded_hubbard8/`. All30 preceding two-spectator obstructions
are integrated into ENERGY v15 and FAMILY v11. Both lower bounds strictly
exceed the freshly replayed preceding pair-family ceilings:
- W0 periodic lower -.6429438983437589; new family cap -.6429248891913111;
  gap 1.9009152447823635e-5;119 positive sources; old-family excess
  3.323439306711367e-5. Open lower -.6429463983437589; physical upper
  -.6106763470511881.
- W1 polished periodic lower -.6607111946158047; new cap -.6551630648953223;
  gap .005548129720482485;118 positive sources; old-family excess
  .00013728253558878796. Open lower -.6607156946158047; physical upper
  -.6184244823693281. W1 numerical limit remains UNRESOLVED.
No claim of exact family attainment. Family caps are not physical energy uppers.

Selected W0 lower comes from W_zero/ and cap from W_zero/priced/.
Selected W1 lower comes from W_plus_1/polished/ and cap from W_plus_1/resumed/.
Final directories consolidate them and have fresh accepted energy/family
replays. Root previous_family/ dirs contain fresh preceding v14/v10 replays.
Final independent closures all pass:30 two-spectator,4pair,14one-spectator,
4spin and coherent hopping moments exactlyzero. Averaged signed-charge laws
have exact CLASSICAL Markov extensions; no full quantum extension follows.
Do not claim the prior30 violations persist in the new final mixtures.

Frozen target U5,t1,V1/4,W-1/5: new open lower -.5226432283437589;
physical upper -.4885616802989547. Without only30 new terms lower
-.5230567883437589, exact loss10339/25000000=.00041356. Previous frozen
pair recipe lower -.5230701060554458; improvement .00042687771168696.
All3 recipes freshly accepted under held_out/, without_two/, previous_recipe/;
frozen_transfer_comparison.json accepted. Fixed-recipe comparison, not a cap
for a reoptimized older family. Same geometry/filling/short range only.

Production: new experiments/marginal_two_spectator_hopping.py implements30
canonical labels i,j,k,l,p,q. Exact projected action enters before PSD test.
ENERGY v15 refuses bad labels and old-version use. FAMILY v11 adds explicit
mode requiring full pair hierarchy; all30 moments exactlyzero; fixed new
coefficient fields forbidden. Nine-shape source cap119 (old89 unchanged in
old mode), local94blocks/4096states/maxPSD200 unchanged. Old family caps
cannot accompany v15 energies. Sources preserved in two_spectator/source_before.

Discovery sources: two_spectator_thermal_numeric.py (118coefficients,500spectral
evaluation bound per run), two_spectator_family_numeric.py (119rows,40pricing
rounds,2eigenvectors/block), two_spectator_family_resume.py, frozen_target.py,
target_replay.py, frozen_comparison.py, separation.py and overlap_replay.py
with the common two_spectator_ prefix. All under discovery/. Proposals remain
untrusted until exact replay. W0 used469 thermal evaluations; W1 initial and
polished each500; iteration/budget stops, no numeric optimality assertion.

Initial sampled family fits had trivial cap0. W0 global40round pricing gave
9276candidates and selected119-source cap. W1 initial40round pricing gave
10506candidates and cap-.653901827567574. W1 clipped radius.02 pricing ran
40rounds but stayed0; polished-point sampling also stayed0. Resuming the
10506candidate ledger for40 more global rounds gave14630candidates and the
selected118-source cap. Last reduced eigenvalues stillnegative; no convergence.
The resume helper INPUT cap is13000, OUTPUT cap20000. The14630 output cannot
be blindly resumed. Next attack needs justified pruning preserving accepted
physical columns or improved constrained/proximal pricing. Do not blindly
increase budgets and repeat unstable global pricing. This is an unresolved
algorithmic choice, not a proved obstruction. No further consistency probe
on these new mixtures has been implemented; possible new conditions untested.

Validation: full945 tests +102subtests pass in849.76s, one existing calibration
return-value warning. Full suite includes new hopping and matching tests.
Focused44 plus matching1 pass. No production source changes after full suite.
All logs and JUnit in two_spectator/. Collector two_spectator_collect.py
audits26 authoritative current receipts and876 source-hash entries;
130historical receipts:4031 unchanged hash entries and190 preserved snapshot
entries. Construction provenance separately checks1302 source-hash entries;
381 files recorded, including unsuccessful artifacts. Old historical sources
resolved via pair_transfer/source_before, pair_transfer/family_source_before,
and two_spectator/source_before. Historical receipts are not fresh current-code
replays. Central marginal_final_validation.json has two_spectator_energy entry.

General quantum representability, generic molecular/long-range/higher-dimensional
transfer and requested-accuracy scalability remain unproved. Preserve refusal
gates and report the W1 family gap honestly. Continue locally; no cloud setup.

## Latest accepted checkpoint — pair-family caps and thirty new violations

Read `research/marginal_pair_family_limit.md` and
`results/marginal_graded_hubbard8/pair_transfer/family_summary.json`.
Previous goal turn made VERIFIED PROGRESS. Goal ACTIVE, not complete/blocked.
All numerical, replay and test processes are terminal. No GPU or paid resources.
Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

LATEST W0/W1 directories: `pair_transfer/<case>/family/` under the result root.
Lower certificates are unchanged from the preceding selected pair-transfer
energy proofs and freshly replayed in these directories. New FAMILY v10 caps:
- W0 periodic lower -.6429787160554459; ceiling -.642977132736826;
  gap 1.5833186198463312e-6 per site; 81 positive physical sources.
- W1 periodic lower -.6608504916575652; ceiling -.6608484771513935;
  gap 2.014506171612037e-6 per site; 81 positive physical sources.
- Open lower W0 -.6429812160554459, W1 -.6608549916575652; physical uppers
  unchanged -.6106763470511881 and -.6184244823693281, freshly checked.
The fixed enlarged family optimum is bracketed at this finite precision.
No exact attainment or physical ground-energy-upper interpretation of caps.

Production: marginal_diagonal_family_limit.py adds explicit pair_transfer mode,
requiring full spectator hierarchy and adding four source slots (89 vs85 for
nine shapes). It checks all4 pair moments exactly. Old modes/caps/refusals stay.
marginal_range_two_family_limit.py adds FAMILY v10. Matching driver adds
joint_pair_transfer_family_proposal_v1 and still refuses v14 energies with old
family caps. ENERGY v14 unchanged. Fixed pair coefficient fields remain invalid
in family certificates, because constraints are moments. No old gate weakened.

Discovery: pair_transfer_thermal_atoms.py exports LAST unrounded optimizer point
from previous W0 initial/W1 polished histories. Largest Gibbs residuals .01465/
.00636 (including positive fidelity residuals), so never accepted as certificates.
Eight/nine dominant atoms plus deterministic rounding perturbations and4096
anchors gave4932/4980candidates. pair_transfer_family_numeric.py retains87
moment equalities and2fidelity inequalities, zero pricing rounds. Selected81
columns solve all89rows exactly with nonnegative weights; both independent
production replays accept. No residual completion or negative clipping.
Numeric discovery solver row cap89; fraction_free_completion.py stays85 unchanged.
Candidate ledgers and diagnostics retained in family/. About35s per numerical
run, ordinary CPU observations only, not a controlled performance benchmark.

Independent closures freshly accepted:4pair,14one-spectator,4spin and hopping
moments exactlyzero. Averaged signed-charge law has exact classical Markov
extension. No quantum density-matrix extension follows.

NEW obstruction: two_spectator_overlap_probe.py constructs all30 PH-even,
reflection-odd five-site charge-hopping products nonconstant in BOTH spectators.
Labels(i,j,k,l,p,q), i<j, k<l distinct spectators; C=q_k^p q_l^q B(i,j),
Y=C-C_reflection, T=Y_left-Y_right. Odd-distance hopping uses powers(1,1)/(2,2),
even-distance uses(1,2)/(2,1). 60primitive products form30reflection pairs.
Full-Fock CAR/Hermiticity/reflection/PH/spin-flip and periodic cancellation
checked. Independent direct bit swaps tested all30actions on4096states.
Two-site hopping row-sum bound2; disjoint charge-power norms<=1 imply ||T||<=8.

ALL30new moments nonzero in BOTH accepted81-source mixtures. Exact restricted
integer Gram rank: old17hopping directions (3unconditioned including2nearest
profiles +14one-spectator) -> combined47. Restriction uses two-bit off-diagonal
entries in N<=4; old projectors onlyN5/6/7, diagonal/spin/pair operators have no
entries there. Thus all30directions independent modulo the current full family.
- W0 largest abs moment .000547192605072811, label(0,3,2,4,2,2).
- W1 largest abs moment .0014597247697291343, label(0,2,3,4,2,1).
Exact rational moments and47x47integer Gram in two_spectator_overlap.json.
Violations survive symmetry averaging and rule out stationary extensions of
THESE averaged mixtures. They do not prove every optimum violates or guarantee
strict energy improvement from adding the new class.

Next substantive work: integrate/test this30-direction class in energy proofs,
then seek compact matching family caps and re-evaluate transfer. This suggests
ENERGY v15 with118numerical parameters, FAMILY v11 with119scalar rows/sourcecap
(89+30); these are prospective counts, not implemented/accepted yet. Preserve
all coefficient, source, rational-weight and full-Fock validation gates. Do not
claim two-spectator-corrected energy: none constructed this turn.

Transfer U5,V1/4,W-1/5 remains the previous frozen held_out/ experiment; no new
coupling or geometry transfer this turn. Generic molecular, long-range/higher-D
transfer, general representability and requested-accuracy scaling still unproved.

Validation: full914tests+102subtests passed616.22s, terminal exit0; one unchanged
style warning in test_v4_calibration.py::test_world returning a value. Focused
family/version/refusal30tests passed105.18s. Two subsequent two-spectator tests
passed5.47s separately after full-suite collection. Sixteen currentreceipts/
508sourcehashentries verified. family_provenance.json hashes249files.

Historical audit:114previous receipts (97thermal-era +17pair-energy) have3564
unchanged hash references and149 matched to preserved pre-change snapshots.
New snapshots in family_source_before/ (3files) supplement source_before/.
Those historical receipts are not described as fresh replays under currentcode.
Central latest section: pair_transfer_family_limit. Prior checkpoints below are
historical and superseded for W0/W1 selections.

## Latest accepted checkpoint — pair-transfer energy integration and frozen transfer

Read `research/marginal_pair_transfer_energy.md` and
`results/marginal_graded_hubbard8/pair_transfer/summary.json`.
Previous goal turn made VERIFIED PROGRESS. Goal ACTIVE; not complete/blocked.
All local numerical/replay/test handles are terminal. No GPU or paid resources.
Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

SELECTED final paths under `results/marginal_graded_hubbard8/pair_transfer/`:
- W0: W_zero/ (initial search; polished rounded proposal was worse).
- W1: W_plus_1/polished/ (improved and independently accepted).
- New coupling transfer: held_out/, plus held_out/without_pair/ ablation.

Exact accepted results (per site):
- W0 periodic lower -.6429787160554459; open lower -.6429812160554459;
  physical upper -.6106763470511881. Strict excess above the entire old fixed
  spectator-family ceiling = 7.922948731852697e-5.
- W1 periodic lower -.6608504916575652; open lower -.6608549916575652;
  physical upper -.6184244823693281. Strict excess above the old ceiling =
  7.080648305529027e-5.
- Held-out U5,t1,V1/4,W-1/5: open lower -.5230701060554458, physical upper
  -.4885616802989547. Projectors, penalties, corrections and filter frozen from
  INITIAL W0 pair certificate; only physical profiles and scalar threshold changed.
  Removing only pair terms gives open lower -.5473074060554458. Exact fixed-
  coefficient improvement .0242373. No ceiling on reoptimized no-pair family.
  Same chain geometry and half filling; no generic molecular transfer claim.

Production changes: experiments/marginal_pair_transfer.py has4 reflected
five-site pair-transfer directions (pairs01,02,03,12), translated difference
on6sites. Independent all-Fock bit-swap tests check signs, PH/spin symmetries,
full-image projection and cancellation. marginal_projector_extendibility.py
ENERGY v14 includes pair matrices before exact PSD, preserving all v13 gates.
Old versions explicitly refuse pair fields. Local coverage unchanged94blocks,
4096states,maxPSD200; opening cost unchanged. Full independent replays accepted.

FAMILY v9 is NOT extended. marginal_diagonal_family_limit.py refuses fixed pair
coefficients; range_two_family_limit_replay.py refuses pair-energy/v14 matching.
The old caps were freshly re-evaluated in each selected previous_family/ using
current code. strict_family_separation.json matches the fixed target, size,
projector sources, ratio, ceilings and sparse span, then proves new lower exceeds
old cap. It does not claim a cap or numerical optimum for the enlarged family.
Fresh pair_transfer_overlap.json in selected previous_family/ rechecks allfour
nonzero moments, exactrank4 and norm<=4 in the OLD ceiling mixtures.

Discovery pair_transfer_thermal_numeric.py extends the old fullmetric softmin
chart from84to88coordinates. Existing signed±8 bounds, other±2,2temperatures,
120iterations/stage,<=500full-spectrum calls. Initial calls334W0/451W1;
polish301W0/412W1. Initial stages and final polish status hititerationlimits.
No optimizer flag accepted. W0polished lower -.6429789239328517 worse; retain
as nonaccepting failed improvement. W1initial lower -.6608529108357133 accepted
but superseded by polished. Both current directories retain dual_atoms.json
(untrusted), thermal_history.json and proposal provenance for next construction.

Next substantive work: construct/certify a compact upper limit for the enlarged
pair-transfer family, then identify any remaining quantum consistency failures.
Current FAMILY v9 has85rows/sourcecap85; adding4pairmoments suggests89rows/cap89,
which requires explicit production extension and tests, not relabeling oldcaps.
Existing fraction_free_completion.py row bound85 also needs deliberate handling
if reused on89rows. No new pair-family mixture is accepted yet. The current
pair-family numerical optimum remains unproved despite strict energy improvement.
General representability and requested-accuracy scalability remain unproved.

Validation: full suite898tests+102subtests passed620.84s, terminal exit0.
One pre-existing-style warning in test_v4_calibration.py::test_world returning a
value; no failure. Focused pair/spectator/spin/hopping67tests passed177.59s.
Additional9separation+oldfamilydriver-refusal tests passed9.95s after full-suite
collection. Seventeen current receipts /613sourcehashentries verified;151files
in provenance manifest. All selected lower/physical upper and oldfamily caps
freshly re-evaluated; exact24site vs160bit upper enclosure included.

Source snapshots and historical evidence: production/source edits invalidate
CURRENT hashes of some earlier receipts. Pre-change4files preserved under
pair_transfer/source_before/ with original_sha256.json. Previous97receipts/
3100hashentries audited:2977unchanged,123match preserved original source. They
remain historical evidence; NOT fresh current-code replays. Current acceptance
is the17explicit newreceipts. Central latest section is pair_transfer_energy;
older sections/checkpoints below remain historical, not current selections.

## Latest accepted checkpoint — thermal limits certified below 2e-7/site

Read `research/marginal_spectator_thermal_limit.md` and
`results/marginal_graded_hubbard8/spectator_hopping/thermal_limit_summary.json`.
Goal ACTIVE. Numerical-limit question resolved to stated finite precision for
fixed W0/W1 family; exact attainment and broader goal remain unresolved.
No GPU/paid resource used. Never contact archived task
01a08e58-e30c-7800-8df6-ab713cbd7787.

LATEST W0/W1: `spectator_hopping/<case>/thermal/final/` under result root.
W+/-0.1 still use polished/. New independent production ENERGY v13/FAMILY v9
replays accepted; no production change or weakened gate. Both77sources.
- W0 periodic lower -.6430581318116138, family ceiling -.6430579455427644,
  gap1.8626884943302653e-7, gap reduction307.8407774x vs completed_basis/.
- W1 periodic lower -.6609213909597719, family ceiling -.6609212981406204,
  gap9.281915150973713e-8, gap reduction4999.5408094x.
- Million-site OPEN lower W0 -.6430606318116138; W1 -.6609258909597719.
  Fresh physical uppers unchanged W0 -.6106763470511881; W1 -.6184244823693281.
Family ceilings bound attainable lower certificates, NOT physical energy.

Discovery thermal stages: spectator_thermal_numeric/polish/scaled/expanded/
fullmetric.py, each<=500full-spectrum evaluations; signed search bounds±2→±8
only in expanded/fullmetric, no production change. Full-spectrum softmin,
Gibbs gradients, diagonal/full spectral Hessian preconditioning. Some stages
hit budgets; success flags NOT certificates. Direct Gibbs residuals nonzero.
Selected lower W0fromthermal/scaled; W1fromthermal/fullmetric. Rounded W0
later trials were worse. All candidate histories retained, nonaccepting.

Family atoms from LAST unrounded fullmetric point:8dominant eigenvectors.
Fullmetric atom maximum residualsW0~4.50e-4,W1~4.79e-5. Rounded perturbations
plus previous sources and4096determinants make4900columns, zero pricingrounds.
Existing spectator_atom_family_numeric.py succeedsW1 with77exact sources.
W0original77basisfailsall85equations. Failed finaldiagnostic retained.
W0successful proposal comes from thermal/affine_completion/, copiedtofinal.
Ledger thermal/ledger/candidate_ledger.json.gz; oldpool completed_basis/85sources.
Eight residual equations,16bounded residualLPs, selectedtwo dependent oldcolumns
rank1. Strict residual_basis_completion refuses; affine_completion retainsone
freeparameter and exactnonnegative interval endpoint. All16exactfeasible;
selectedoldindices45,62; zero oldselectedcoefficients removed yields77sources.
All85originalequations/trace/nonnegativity/sourcecap/weightlimits rechecked.
New discovery spectator_residual_affine_completion.py preserves failedvariant.

Fresh pair tests allfour violated, exactrank4/norm<=4/CAR/symmetry/cancellation:
W0moments pairs01,02,03,12 = +.0011728671459324202,-.00016328053176235783,
+.00028588480533632224,+.000014527907031592661.
W1 = +.00040968890751056253,-.00013165713093398864,
-.0009752628901861763,+.000030488303052099796.
All14spectator/4spin/unconditionedhopping zero; averaged signed-charge law has
classicalMarkovextension. No quantumextensionclaim. Pair-transfer energy
integration remains next unresolved experiment; violation of THESE mixtures
does not prove every optimizer violates or guarantee strictenergyimprovement.
No fresh ablation atnewlower coefficients. No general molecular/longrange/
higherD transfer or requested-accuracy scaling theorem.

New spectator_interval_refinement.py checks samefixedfamily, exact nestedbounds,
strictgapreduction, inputreceiptprovenance.25focusedtests passed0.07s (16new
intervalcomparison+9unchangedfractionfree); no fullsuitererun. Previous847tests
+102subtests historical only. Sixteennewacceptingreceipts,97total audited,
3100sourcehashentries current. thermal_limit_provenance.json hashes370files
includingallrelevant imported discoverymodules/physicalmodules andfailedtrials.
Collector spectator_thermal_collect.py updates central spectator_thermal_limit.
No controlledperformanceclaim; ordinary CPU observations only.

Previous checkpoint below is historical and superseded for W0/W1.

## Latest accepted checkpoint — exact completion recovers both rejected ceilings

Read `research/marginal_spectator_basis_completion.md` and
`spectator_hopping/basis_completion_summary.json` under the result root.
Previous goal turn made verified progress. Goal ACTIVE; not complete/blocked.
No running local calculation or GPU. No GPU use/spend this turn. Never contact
archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

LATEST W0/W1 paths: `spectator_hopping/<case>/completed_basis/`.
W+/-0.1 still use prior polished/. All lower/physical upper endpoints unchanged
from limit_refined; both were freshly replayed at the new paths. New ceilings:
- W_zero: periodic ceiling -0.6430028839174731; remaining gap 5.7341147414e-05; ceiling improvement 1.0566025428e-05;85sources.
- W_plus_1: periodic ceiling -0.6605408228122668; remaining gap 0.000464053135867; ceiling improvement 4.70179513636e-05;85sources.

The exact-basis blocker for the previous numerical proposals is recovered.
Do NOT equate this with the unrestricted family optimum: gaps remain5.73e-5
and4.64e-4, and pair-transfer obstructions persist. New W0 witness uses85sources
vs prior77; W1 stays85. No productionverifier or source/parsercap changed.

New discovery construction:
- fraction_free_completion.py transforms integer selected columns A, target b,
  and old completion pool B to D*x+P*z=r, C*z=s. All divisions exact, row swaps
  checked; original equations are verified again after reconstruction.
- Raw source columns use integer amplitudes and norm multiplication to clear
  source denominators. LCM scaling of fidelity rows includes projector norms,
  ratio denominator and both ceiling denominators. All85 equations retained.
- Current LP basis from saved ledgers:83columns W0,84W1. Completion pool is the
  OLD polished85-source strictlypositive mixture, not the latest77source W0.
- Enumerate all3570old-column pairs forW0,85single columns forW1. Nonnegative
  completions14/1. W0selectedoldindices56,81; added normalizedweights about
  1.3746911336272823e-12,3.0688246411071267e-12. W1selectedoldindex39, determinant
  {'21':1}, weight2.804448662337289e-10. Exactfinal85source constraints accepted.
- Original spectator_basis_completion.py succeedsW0. W1oldpool includes two
  individualdeterminants (indices39/40) outside a single invariantblock; first
  run refused. completion_before_anchor.log retained. New
  spectator_basis_completion_anchors.py evaluates determinants directly,
  including spin diagonals, verifiesoldmixturemoments, then succeedsW1.
- These runs use scaled_basis/candidate_ledger.json.gz (8346/8820states).
  No spectral sampling/pricing rerun. Source/fixedfamily/ledger matching checked.
- New ceilings are about1.63e-11/3.79e-10 ABOVE old numericalLPvalues.
  W0completion can add vectors outside originalledger: no finite-original-LP
  optimality claim. All proposaldiagnostics remain accepted:false untilreplay.
- Bounds:85rows/85completioncolumns/4000combinations/input4096bits/
  pivot50000bits. Actualpivots4060/4284bits; maxweightstringlength2429/2497
  below4096parserlimit. Reduction14/18s, total31/34s; concurrentobservations,
  not controlledperformance/scalability comparisons.

New final accepting receipts percase (8each): range_two_replay,
range_two_family_limit_replay,ceiling_improvement,pair_transfer_overlap,
one_spectator_overlap,spin_overlap,coherent_overlap,charge_markov_extension.
ceiling_improvement checksidenticalfixedfamily andsameperiodiclower thenexact
strictpositiveceilingdifference. All14spectator/4spin/unconditionedhop moments
zero; averagedcharge law hasclassicalMarkovextension. No quantumextensionclaim.

Pairmoments pairs01,02,03,12:
W0 +.0012010158000890039,-.0002063094775837515,+.00023619857827046254,
-.000026958705637268918.
W1 +.0010077858534148759,+.00019815252240568644,-.001125311395993989,
-.00021383178188380064.
All four nonzero, independent exactrank4/fullCAR/symmetry/cancellation/norm<=4
checks rerun. Pairtransferenergy integration still not done.

Validation:9focusedtests pass0.03s in tests/test_marginal_fraction_free_completion.py.
They cover originalequations,rowswaps,negativedeterminant,inconsistentrectangular
systems,residualcompletion,anddependent/noninteger/oversizedrefusals. Lasttestlog
basis_completion_focused_validation.log; hashes in focusedvalidationJSON.
16new acceptingreceipts,81totalcurrent/historical,2650sourcehashentries match.
Prior847test+102subtest fullsuite NOT rerun; production ENERGYv13/FAMILYv9
unchanged. Centralstatus section spectator_basis_completion. Provenance and
collector: basis_completion_provenance.json / spectator_completion_checkpoint.py.

Next: use this completion method with better numericalcandidate sets, work on
remainingprimal/dual gap, or integrate4verifiedpairtransferdirections. One
untested numericalidea is global spectral soft-min smoothing/annealing with
full-Fock eigenvalue checks, since prior SLSQP is nonsmooth and boundedconic
models can return invalidobjectivevalues atuser_limit. Do not treat thatidea
as executed orverified. Existing ledgers remain reusable. General quantum
representability, molecular/longrange/higherD transfer and requested-accuracy
scalability remainunproved. Preservefullgoal andleaveactive.

## Latest accepted checkpoint — tighter W=0 limit, numerical optimality still open

Read `research/marginal_spectator_limit_refinement.md` and
`spectator_hopping/limit_refinement_summary.json` under the result root.
Previous turn made verified progress. Goal ACTIVE, not complete or blocked.
No live calculation or GPU remains. No GPU use/spend this refinement. Never
contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

NEW selected W0/W1 artifacts: `spectator_hopping/<case>/limit_refined/`.
W+/-0.1 remain at their prior `polished/` checkpoints. Historical spectator
combined_summary.json retains old values; latest overrides are separate.
- W_zero: exact open lower -16076568126622177/25000000000000000 = -0.643062725064887; upper -0.610676347051188; gain 1.46592110928e-06; periodic family gap 6.79071728421e-05; sources 77.
- W_plus_1: exact open lower -826261719935167/1250000000000000 = -0.661009375948134; upper -0.618424482369328; gain 1.0178688312e-06; periodic family gap 0.00051107108723; sources 85.

W0 family gap is ~3.09x smaller with77sources instead of85. New physical uppers
freshly computed, unchanged. Production ENERGYv13/FAMILYv9 unchanged. All94
localblocks/sum4096/max200/151PSDchecks and all prior constraints remain.
Selected lower origin for BOTH targets: subspace48 (12 rounds, maxcone48,
8activeblocks, temporary radius.05, nominal15solver seconds/60iterations each).
Neither converged; exact full-Fock replay accepts their rounded best candidates.
Selected W0 family origin: subspace48/diagonal_family_limit_proposal.json from
spectator_atom_family_numeric.py, scale.0002, 200extra vectors from8dominant
conicdual eigenvectors, determinant anchors, ZERO pricing rounds.77sources.
Selected W1 family origin: prior polished85source proposal, freshly matched.

Methods and failure evidence:
- spectator_cutting_numeric.py: supporting-plane LPs, tiny exactacceptedgains;
  both ended HiGHSUnknown, NOT numericaloptimum.68/334matrixevaluations.
- subspace24 W0 stoppeddimensioncap withoutimprovement. Initialsubspace48
  tuplekeyJSON exportfailed; fixed, preservedpre-fixsource/logs, rerunsaccepted.
  Two dependentatomexportcalls failedmissingfile, then reran aftersourceexisted.
- Frozen temporary-box dual mixtures have nonzero full-family residuals:
  max~1.85e-5 W0,2.24e-4 W1. Eightdominantdualatoms were only candidatevectors.
  Their tight-looking model values do NOT certify unrestrictedfamilyoptimum.
- global_subspace removed temporarybox, initializedall94cones; stoppedtotal
  compressedorder1000, noimprovement. global_diagonal additionally enforced
  all4096 projected diagonal inequalities, cap2000/maxcone48/12rounds.
  No selectedimprovement; finaluser_limit values can lie BELOW acceptedlower.
  New stdlib spectator_solver_limit_refutation.py proves exactcontradictions
  for2W0 and3W1reported decimalvalues. These values cannot be familyceilings;
  not a refutationofduality. Actuallastsolvertime~16.29/15.46s vsnominal15.
- atom_pricing 40rounds/twoeigenvectors: W0 basis83, W1basis84, exactinconsistent.
  Numericalcaps ~-.6430028839338047 /-.6605408231910822 NOT accepted.
- scaled_basis: finalLP row scales100,10000,1000000, originalexactequations;
  all exactbasesinconsistent. Both finitecandidate sets saved in
  scaled_basis/candidate_ledger.json.gz. Reuse these instead of redoingpricing.
- slack_repair readstheledgers and triesadding omittedfidelity slack columns.
  W0bothslacks: dependent; single slack: negativeexactweights; noneinconsistent.
  W1single slack: negativeweights; noneinconsistent. No proposalaccepted.
  Furtherrepair needsa DIFFERENT basis, not weakenedmomentchecks/clippedweights.

New final receipts (10each): range_two_replay,range_two_family_limit_replay,
family_escape,spectator_ablation,pair_transfer_overlap,one_spectator_overlap,
spin_overlap,coherent_overlap,charge_markov_extension,solver_limit_refutation.
All exactaccepted. Old14spectator/4spin/unconditionedhop moments zero.
Charge projection exact stationaryMarkovextension remains classicalonly.
All FOUR pairtransferobstructions persist; next energyintegration remainsundone.
W0 moments pairs01,02,03,12:
+.001194488138198365,-.00020492360144159502,+.0002564377626730751,
-.000026976852488055283. W1 prior85mixture unchanged (moments inoldcheckpoint).
Pairrank4/fullCAR/symmetries/periodiccancellation/norm<=4 rerun.
All oldspin-family strictescapes and currentfixedcoefficient ablations rerun.

Validation:20newfinal +5intermediate +40prior =65receipts;2200sourcehashentries
match. No productionchange, prior847tests+102subtests regression NOT rerun.
Fresh exactreplays andindependentphysicalprobes are thisrefinement's evidence.
Collector: discovery/spectator_limit_checkpoint.py; centralstatus section
spectator_limit_refinement; provenance manifest limit_refinement_provenance.json.
All executionhandles terminal andfinalprocessscan foundnone.

Next options: exactbasis recovery using savedfinitecandidateledgers; better
bounded conic optimization preservingfullphysicalgates; or integrate4verified
pairtransferdirections. Do not pretendcurrentfamilyoptimum resolved. General
representability, arbitrary molecular/longrange/higherD transfer and requested-
accuracy scalability remainunproved. Preservefullgoal andleaveactive.

## Latest accepted checkpoint — spectator hopping escape, pair-transfer obstruction

Read `research/marginal_spectator_hopping.md` and
`results/marginal_graded_hubbard8/spectator_hopping/combined_summary.json`.
This supersedes the spin checkpoint below. Goal ACTIVE; verified progress,
not completion. No local calculation or GPU remains running. No GPU use/spend
this milestone. Never contact archived task 01a08e58-e30c-7800-8df6-ab713cbd7787.

Final artifacts: `spectator_hopping/<case>/polished/`. All four million-site
open-chain U4,t1,V1/2 exact intervals and matched family ceilings accepted:
- W_zero: lower -16076604774649909/25000000000000000 = -0.643064190985996; upper -0.610676347051188; strict periodic escape 6.23508252834e-05; current periodic family gap 0.00020978666201.
- W_plus_1_10: lower -16092802273917671/25000000000000000 = -0.643712090956707; upper -0.611451160583002; strict periodic escape 5.68915972403e-05; current periodic family gap 0.00022013538764.
- W_minus_1_10: lower -16066401424763747/25000000000000000 = -0.64265605699055; upper -0.609901533519374; strict periodic escape 9.8500088993e-05; current periodic family gap 0.00014669289457.
- W_plus_1: lower -413131496135603/625000000000000 = -0.661010393816965; upper -0.618424482369328; strict periodic escape 8.58084851227e-05; current periodic family gap 0.000512088956062.

Strong conclusion: the new periodic lower strictly exceeds the ENTIRE previous
fixed spin+hopping+complete-charge family ceiling for every target, keeping
sources, ratio, ceilings, nine-shape sparse span and physical target matched.
Each family_escape.json freshly replays the old spin witness with current
production code. This does not cover other sources or arbitrary certificates.
Physical uppers freshly recomputed, exact 24-site contraction inside 160-bit
interval enclosure, million-site recurrence; unchanged from prior milestone.

New independent class: B(i,j)=sum_spin(c_i†c_j+c_j†c_i), C=q_k^p B with k distinct,
p=2 for odd distance and p=1 for even. Five-site Y=C(i,j,k)-C(4-j,4-i,4-k),
six-site T=Y_left-Y_right. Fourteen canonical triples: nine q², five q.
Complete nonconstant one-spectator charge-function real-hopping class after
PH-even/reflection-odd projection. Translation cancellation adds no physical
coupling/opening penalty. Gram rank17 combines old3 unconditioned directions
with new14 on N<=2; prior diagonal/spin entries absent there and actual fixed
projectors have N5,6,7. Parent basis_probe.json files retain OLD spin witnesses:
all14 moments nonzero. Do not overwrite them with current closure receipts.

ENERGYv13 adds bounded exact spectator_hopping dict; old versions reject field.
FAMILYv9 adds14 exact zero moments, source cap85 (actual85 every final witness).
Old caps and every previous gate retained. Input refuses fixed spectator
coefficients; matching driver requires v9. New helper:
experiments/marginal_spectator_hopping.py. Independent CAR oracle/probes in
results/marginal_graded_hubbard8/discovery/spectator_hopping_basis.py.
94 local blocks / sum4096 / maxPSD200 /151totalPSDchecks /64sparsecap unchanged.

Numerical primal chart85: p,q,alpha,beta,9sparse,52signed (0:65), gamma65,
spin66:70, spectator70:84, ell84. Coherent array19=[hop,4spin,14spectator].
All first AND polish passes hit240iterations (240/241matrix calls; hard250cap).
No convergence claim. Fresh affine discrepancy <=3.56e-15. Current enlarged
family gaps remain LOOSE; especially W1 ~5.12e-4. Optimization is still open.

Dual85 rows: trace0, nearesthop1:3,61diagonal3:64,19coherent64:83,fidelity83:85.
Numerical directions82. Determinant anchors include spin diagonals.
Exact rectangular solve retains every equation; LPtol1e-10/support1e-14.
Candidate caps1200sample+4096anchor,13000total. One-eigenvector <=80rounds;
two-eigenvector <=40rounds preserves total cap. Final provenance manifest:
- W0 / W+.1: parent proposals, 1eigenvector60rounds, 8130/8008 candidates,
  copied into polished and freshly matched/replayed; preserved driver is
  pre_spectator_sources/spectator_family_before_multieigen.py.
- W-.1 / W1: polished seed, current driver, 2eigenvectors40rounds,
  8482/8470 candidates. W-.1 recovered from initial exact84-column
  inconsistency; W1 improves ceiling. Seed and rounds changed too, so not a
  controlled eigenvector-count comparison. No conic/GPU rerun.
- Initial four ablations launched before required energy receipts existed and
  failed FileNotFoundError. Logs preserved as ablation_before_energy_ready.log;
  all four rerun successfully after accepted energy. No pending failed job.

All polished receipts: range_two_replay, range_two_family_limit_replay,
family_escape, spectator_ablation, pair_transfer_overlap, one_spectator_overlap,
spin_overlap, coherent_overlap, charge_markov_extension (.json each).
Current14 spectator moments exactlyzero by independent CAR; old4spin and
unconditionedhop close. Removing ONLY spectator causes exactnegative residual
in ablation, same physical integer vector nonnegative when restored.
Averaged charge law has exact order5 stationary classical Markov extension,
all243prefix/suffix equalities and729 reconstructions. Not quantum or globalN.

NEXT EXACT OBSTRUCTION: pair transfer, not yet energy-integrated.
d_i†=c_i_up†c_i_down†, Jij=d_i†d_j+d_j†d_i. Four canonicalpairs01,02,03,12
with reflected five-site differences / six-site translates. All four moments
nonzero in each current family, PH/reflection invariant. NormJ=1 via independent
16state check, triangle telescope norm<=4. Rank4 on two-electron double states:
[[8,0,0,-4],[0,8,0,0],[0,0,12,0],[-4,0,0,12]]. Independent modulo current
family: pair moves4bits vs singlehop2; spin cannot move double; projectors vanish.
Strongest moments W0 +.0012134183152451856 (01), W+.1 +.0008866224579082949
(01), W-.1 +.0015644050340553905 (01), W1 -.0012895082622640165 (03).
Exclude stationary quantum extension of averaged six-site density matrix;
not all inhomogeneous extensions or general representability. Remaining
multi-spectator/end-point-conditioned/spin correlations outside current class.

Validation: focused28passed111.41s; FULL847passed+102subtests661.31s, one
existing tuple-return warning. JUnit full_validation.xml:949 including subtests,
847 testcase elements,0errors/failures/skips. Full run completed exit0 before
report/checkpoint collector added; collector subsequently executed successfully.
40 mathreceipts,1348sourcehash entries allmatch. Latest central status updated.
No extra suite needed for report-only changes; no production changes since run.

Next work must address the loose numerical limit and/or integrate the four
verified pair directions into compact exact certificates. General quantum
representability, molecular/long-range/higher-D transfer and requested-accuracy
scalability remain unproved. Do not equate a finite family ceiling with a
physical energy upper or mark the broader goal complete.

## Latest accepted checkpoint — strict spin-family escape, conditioned hopping obstruction

Read `research/marginal_spin_telescope.md` and
`results/marginal_graded_hubbard8/spin_telescope/combined_summary.json`.
Supersedes the hopping checkpoint below. Goal ACTIVE; this turn made verified
progress. No local jobs or GPU remain; no GPU ran/spend this turn. Never contact
the archived old task.

BEST final artifacts are `spin_telescope/<case>/polished/`.
All4 exact million-site OPEN U4,t1,V1/2 intervals:
- W_zero: open lower -1607889346956211/2500000000000000 = -0.6431557387824844; upper -0.6106763470511881; gain 2.60054126252e-05; strict old-family escape 6.42480255042e-06; new periodic family gap 2.91969712046e-05.
- W_plus_1_10: open lower -2011821275356161/3125000000000000 = -0.6437828081139715; upper -0.6114511605830021; gain 4.59073600913e-05; strict old-family escape 3.34128650956e-05; new periodic family gap 1.38255600244e-05.
- W_minus_1_10: open lower -4017373365675113/6250000000000000 = -0.6427797385080181; upper -0.6099015335193740; gain 2.3169285742e-05; strict old-family escape 9.63185772429e-06; new periodic family gap 2.51814284752e-05.
- W_plus_1: open lower -16531480413486909/25000000000000000 = -0.6612592165394764; upper -0.6184244823693281; gain 0.000364535844344; strict old-family escape 0.00033350036236; new periodic family gap 0.000163014237389.

All4 first and polished energies accepted, physical uppers freshly recomputed
with exact24-site contraction in160bit enclosure; uppers unchanged. All4
polished family and strict escape receipts accepted. ParentW0 family accepted;
parentW+.1 family numeric failed, so don't use parent dirs as final families.

STRONG result: each new accepted periodic lower STRICTLY EXCEEDS the previous
hopping+complete-charge family ceiling, with the SAME sources/ratio/ceilings/
sparse span/physical target. `family_escape.json` freshly replays the old
physical mixture using current production code, matches family fixed fields
and accepted new energy hashes, and proves positive difference. This goes
beyond all old-family scalar tuning; no claim of universality/changed sources.

Production ENERGYv12: `spin_telescope` exact bounded dict of1..4 canonical
pairs 0,1;0,2;0,3;1,2, nonzero afterpruning. DefineAij=4Si.dot.Sj,
Y=Aij-A_(4-j,4-i) on5sites, T=Yleft-Yright on6. Newhelper
`experiments/marginal_spin_telescope.py` buildsCAR spin-dot actions andchecks
all4096 Hermiticity/spin-number/reflection, thenexistinghop.projected_matrix
projectsE^TTE. v12retainsallpreviousfields andoldversionsrejectspinfield.
SpinT cancels undertranslation: no physicalspin coupling or openingcost.
All94localblocks/maxPSD200/4096total/151totalPSDchecks/64sparsecap unchanged.

Production FAMILYv8: shared private_replay `spin_telescope=True`requires
hoppingmodeandallpreviousflags; adds4sourcecap and4exactzeromoments. Cap71,
actual71everyselectedwitness. Oldcapsunchanged. Allprior52signedpatterns,
profiles, compactquadratic/indicator/sparse/hopping andfidelityconditions
remainexact. Matchingdriverrequiresv8forspinenergy andchecksallsources/shape
span/ceilings/ratio/target. Newcapscoverunrestrictedspincoefficients butNOT
conditionedhopping/newspaces. Newfamilygaps~1e-5 exceptW1~1.63e-4: CURRENT
SPINNUMERICALOPTIMUM stillnotsharplyresolved, evenwhereSLSQPflagsuccess.

`spin_span_rank.py` exact4x4Gram in(1,1)two-electron double-spin-flip entries:
[[32,0,0,-16],[0,32,0,0],[0,0,48,0],[-16,0,0,48]], rank4.
Olddiagonalentriesvanishthere; hoppingchanges2bitsnot4; actualfixedprojectors
verifiedin(3,3),(2,3),(3,2),(3,4),(4,3), vanishin(1,1). Thusfourindependent
newdirections moduloentirepreviousfixed-source certificatespan, penalties
included. Receiptrootspin_rank.jsonincludesactualpolishedW0certsourcehash.

Discovery `spin_numeric.py`: oldsignedprepare physicalbaseEXCLUDESsigned,
hoppingtelescopeandspin; physicalbaseincludesfixedoldquad/square/indicator/
range2 profiles. Blockcoherentarrays have5matrices [hopT,spin01,spin02,spin03,
spin12]. Chart71 vars: p,q,alpha,beta,9sparse,52signed,gamma,4spin,ell.
Indices gamma65,spin66:70,ell70. Firstsearchesallhit240iterations/241matrix
calls. Onepolishpass W0success95calls, W-.1success4calls; W+.1limit240calls,
W1limit241calls. Perpasshard250matrixcalls unchanged. FreshCAR difference
<=3.56e-15. No convergence-rate orperformanceclaim.

`spin_family_numeric.py`:71exactrows trace0,nearesthop1:3,61diagonal3:64,
5coherent64:69,fidelities69:71. Numericaldirections68. IMPORTANT: determinant
anchors'4spin diagonalexpectations areNONZEROin general (hopTdiagzero).
Currentdriver includesexactspinanchor values. LPtol1e-10,support1e-14,
--determinant-anchors --pricing-rounds60; cap80rounds/13000candidates; standard
71-rowrectangularsolve retainsallrationalequations/positivity. W0parent7334,
W-.1parent7540,W1parent7566candidates; theirproposalscopiedtopolishedthenfresh
matchedreplay. W+.1parentfailed70-columnexactbasisconsistency at7426candidates;
sameboundedsettingsonpolishedseed recovered71-columnproposal at7326candidates.
No gate weakened. Parentfailurepreserved; selected_numeric_provenance.json
recordsactualseed/diagnosticpaths,currentdriverhashandselectedproposalhash.

`spin_ablation_probe.py`: all4polishedspin_ablation.json accepted physical
integerRayleighvectors showdeletingONLYspinT atsubmittedcoefficients breaks
localPSD, whilewithspinresidualnonnegative. HopT retained. Acceptancecontracts
freshdeterminantCARactions, notnumerical blockmatrices.

NEW NEXT physical obstruction: `conditioned_hopping_probe.py` tests9canonical
(i,j,k) on5sites, i<j odddistance, kdistinctspectator, tuple<reflectedtuple.
C(i,j,k)=q_k^2 B(i,j), B=sumspin(c†c+h.c.); Y=C-Cref, T=Yleft-Yright.
All9exactmomentsNONZEROinall4spin-completefamilywitnesses. Strongestnormalized
(conservativecommon normupper8): W0 (0,3,1) .0023259915105646144;
W+.1 (0,3,1) .00225723035773537; W-.1 (0,3,1) .002483081628874769;
W1 (1,2,0) -.0035739023645967893. Exactall4096Hermiticity/spin-sector/PH/
reflection/spin-flip invariance andperiodicsumzero checked. p=q² commuting
spectatorprojector norm1, Bnorm<=2, triangleTnorm<=8; NOTexactTnormclaim.
Nonzero momentspersist afterPH/reflectionaveraging andforbid stationaryquantum
extensions matchingaveraged6sitemarginals; notinhomogeneoussingleextension.
No conditioned-hopping energyintegration/optimization orindependencerank yet.
Nextboundedattack: exactrank/representation thencompactconditionedhopenergy+
matchingfamily, or workoncurrentloosespin-family numericalgap. Allnine are
TESTEDdirections, notclaimofcompleteallquantumconstraints. Othercharge/sign-
conditionedhopping andhighercorrelations remainbeyondthisscan.

All4old`spin_overlap.json` now exactzero forfourspinmoments, and
`coherent_overlap.json` exactzero forhopT, independentlyreconstructed.
`charge_markov_extension.json` all4still exactPH/reflectionAVERAGED classical
signed-chargeorder5stationaryextension (243prefix/suffix,stochastic,stationary,
729reconstructedprobabilities). NOTspin/coherence/globalfixedNrepresentability.

VALIDATION:14newspintests. First12energy/helpertests16.45s; focused27spin+hop
passes46.14s after2newfamilytests. Fullrepositoryprocess exited0 andreported
`824 passed, 1 warning, 93 subtests passed in 548.25s`. Subsequentread-only
collectionlists833items. Preservecountsverbatim; do notinfer missingfailures
orclaim833reportedpasses. Onewarninguntouchedtest_v4_calibration::test_world
returnstuple. Logfull_validation.log, collection.log, centralsummaryupdated.
All32selectedacceptingreceipts(8percase)plusrankreceipt have1107source/input
hashentriesverified. Sourceversionsbeforestage savedpre_spin_sources/manifest.
Parentfirstacceptedenergiesretained, but polished paths are currentbest.
No localresearch/testprocess remainsrunning. Generalrepresentability,
arbitrarymolecular/longrange/higherDtransfer,requestedaccuracyscalability
remainunproved; finitechainresultdoesnotsolvegeneralquantumchemistry.

## Latest accepted checkpoint — compact hopping, remaining spin obstruction

Read `research/marginal_hopping_telescope.md` and
`results/marginal_graded_hubbard8/hopping_telescope/combined_summary.json`.
This supersedes the earlier signed-charge checkpoint. Goal ACTIVE; previous
turn made verified progress. No local research/test jobs remain, no GPU ran
or spend occurred. Never contact the archived old task.

Current best million-site OPEN U4,t1,V1/2 targets:
- W_zero: open lower -16079543604877741/25000000000000000 = -0.6431817441951097; upper -0.6106763470511881; gain 3.409954162e-07; periodic family gap 1.95806100748e-05.
- W_plus_1_10: open lower -1609571788685157/2500000000000000 = -0.6438287154740628; upper -0.6114511605830021; gain 6.1155923352e-07; periodic family gap 1.24944949957e-05.
- W_minus_1_10: open lower -8035036347422001/12500000000000000 = -0.6428029077937600; upper -0.6099015335193740; gain 4.1290980604e-07; periodic family gap 1.35374280177e-05.
- W_plus_1: open lower -16540593809595521/25000000000000000 = -0.6616237523838209; upper -0.6184244823693281; gain 1.41963368038e-05; periodic family gap 3.10354819842e-05.

All four lower and freshly recomputed physical upper receipts accepted.
Exact24-site contraction is inside the160bit upper enclosure; uppers unchanged.
Production ENERGYv11 adds scalar `hopping_telescope` gamma (nonzero exact,
old versions reject) for T=B03-2B14+B25, B=sum_spin(c†c+h.c.). It is Y_left-
Y_right with five-site Y=B03-B14. Periodic translated sum zero; no physical
range3 hopping and no extra opening penalty. Helper
`experiments/marginal_hopping_telescope.py` rebuilds12CARwords, checks all4096
states for Hermiticity/reflection/spin-sector invariance, and projects E^TTE.
All94local blocks/max200/4096total/151totalPSDchecks and64sparse cap retained.
Energygamma W0 923/1e6; other exact values in combined summary/certificates.

Production FAMILYv7 adds exact zero hopping moment and one source:67cap,
all67actual for every selected target. All52signed patterns, allprior profiles,
quadratics/indicators/sparse shapes and fidelity inequalities remain exact.
Shared private diagonal replay `hopping_telescope=True` requiresfullsigned;
oldsourcecaps preserved. Matched familydriver refusesv11energy witholdercap.
Current family ceilings cover unrestrictedreal hopping telescope coefficient,
fixedprojectors/ratio/ceilings/nineshapes andallmean-correctprofiles. They do
NOT cover spin corrections or other supports/sources. Exact current gaps
are~1e-5: numerical optimum is still NOT tightly resolved. Three SLSQPflags
success, W1hits240iterations; no numerical success taken asoptimalityproof.
No strict separation from the fully reoptimized older signed-charge family
has been proved; accepted energy improvements also include further tuning.

Discovery `hopping_numeric.py`: imports oldsignedprepare, adds exactprojected
T and gamma coordinate.67variables, hard250matrixcalls; first3converge andW1
240iterations. Freshall94CAR mismatch <=3.56e-15. Exactphysical ablations in
`hopping_ablation.json` prove deletingonlyT at submittedcoefficients makes a
physicalinteger Rayleigh residualnegative; withT residualnonnegative. Probe
uses freshdeterminant CAR contractions foracceptance, no numericalmatrixreuse.

Discovery `hopping_family_numeric.py`: physical low-eigenvector pool plus4096
physicaldeterminantanchors, and bounded LP-directed pricing acrossall94blocks.
--pricing-rounds cap80, selected60; totalcandidatecap13000. Optional radiusclip
around primalmultipliers failedtoimproveW0 andisnotselected. NOpricing optimum
claim evenif clippedrun stops. Initialunrestrictedruns for3nonzeroW failed
exactbasisconsistency (65/66columns). TighteningHiGHStolerance1e-9->1e-10 and
supportthreshold1e-10->1e-14 yielded67columnexactpositiveproposals. Bothsettings
changed; separatecause notisolated. Production neverrelaxesanyrationalgate.
Current code writes `proposal_written` diagnostics; failedruns may leaveprior
proposalpath untouched, so do notpromoteunlesscurrentdiagnosticconfirmsoutput.
Selectednumerical provenance: W0 uses snapshotbefore_stabilization.py; others
currentstrictdriver. `selected_numeric_provenance.json` recordshashes/arguments.
Selected currentcandidates W0 7088, W+.1 6952,W-.1 7040,W1 6998. Sourcevectors
remain <=400support/1e9amp. Rectangularsolve checksall67rows exactly.

NEW remainingobstruction: `spin_overlap_probe.py`, all4`spin_overlap.json`.
Aij=4 Si dot Sj; five-siteY=Aij-A_(4-j,4-i), localT=Yleft-Yright.
Fourcanonicalpairs (0,1),(0,2),(0,3),(1,2), allfour momentsnonzero ineachdual.
Strongestabs moment (conservative commonnormbound12):
 W0 pair03 .020931031603670106; W+.1 pair03 .013968157351026295;
 W-.1 pair03 .02815606312265058; W1 pair01 -.03593976758782963.
All4096CAR Hermiticity,spin-sector,PH/reflection/spin-flipchecks, periodic
translatedcancellation, and independenttwo-site singlet/triplet blockchecked.
Pairnorm3, triangle telescopeuppernorm12; NOTexacttelescope normclaim.
These momentspersist afterPH/reflectionaveraging andruleout stationaryquantum
extensions matchingtheaveraged6sitedensitymatrices. Notinhomogeneoussingle-
marginalextendibility/generalrepresentabilityclaim.

Independent`coherent_overlap.json` all4now exactlyzero forprevioushoppingT.
`charge_markov_extension.json` all4still exactPH/reflectionAVERAGED classical
order5stationaryextension:243prefix/suffix+stochastic+stationary+729reconstructed
probabilitychecks. Doesnotextendspins/coherences/fixedglobalN. W0all243states
and715sixchargepatternspositive; actualcounts foreachincombinedsummary.

NEXTboundedattack: integratecompactfourspin-dot telescopes (or justifiedsubset)
intoenergyreconstruction, withmatchingfamilyzero moments/physicalsourcecap,
focusedacceptance/refusal tests, boundedoptimization thenindependentexactenergy
andfamilyreplay. Currenthopfamilylimits doNOTcovernewspincoefficients. Allfour
spinpairs asrepresentabilityseparators aredone; no spinenergyintegration or
optimization yet. Improvingcurrentloose~1e-5familygap is alsoopen. The previous
conicexit137attempt remainsfailed; no rerun thisturn. Don't blindlyrestart it.

Validation: focused18tests25.82s. FULL REPOSITORY819tests+102subtests491.19s
passed; onewarning fromuntouched test_v4_calibration.py::test_world returning
tuple. Broadercommand `python -m pytest -q` vsprior596marginaltests, so don't
comparecounts asaddedtests. Currentnewfile has13focusedhoppingtests. Log
hopping_telescope/full_validation.log; centralmarginal_final_validation updated.
All24selected acceptingreceipts(energy,family,ablation,hop,spin,Markov x4)768
source/inputhashentriesverified. Changedpre-stage proofsources savedunder
pre_hopping_sources/manifest.json. Unselected/faileddualrun snapshots retained;
canonical selectedreceiptpaths are authoritative. No currentjobs orGPU.
Generalrepresentability, arbitrarymolecular/longrange/higherDtransfer, requested-
accuracyscalability remainunproved; goal staysactive.

## Latest accepted checkpoint — complete signed charges, off-diagonal obstruction

Read research/marginal_full_signed_charge.md. Supersedes older checkpoints.
Goal ACTIVE. Previous turn made verified progress. No local jobs remain.
Never contact archived old task. No GPU ran/spend this turn.

BEST million-site OPEN U4,t1,V1/2 intervals/site:
 W0 [-.6431820851905259,-.6106763470511881], exact lower
 -8039776064881573/12500000000000000; gain4.170439359016e-5.
 W+.1 [-.6438293270332963,-.6114511605830021], exact lower
 -2011966646979051/3125000000000000; gain4.217837253284e-5.
 W-.1 [-.6428033207035662,-.6099015335193740], exact lower
 -16070083017589153/25000000000000000; gain4.112499390368e-5.
 W1 [-.6616379487206246,-.6184244823693281], exact lower
 -3308189743603123/5000000000000000; gain6.563016716760e-5.
All4energyreceiptsaccepted; uppers freshlyrecomputedwith24siteexactenclosure,
unchanged. W0physicalHamiltonianunchanged. MaxlocalPSD200,94blocks,total4096,
151totalPSDchecks andsparse64entrycap unchanged.

NEWrepresentation: marginal_signed_charge_telescope.py PATTERNS52 disjoint
PH-even reflection-odd signed-charge pattern classes on5sites. Keys comma-
separatedqvalues(-1,0,1);eachclass+1onq/-q,-1onreflectedq/-q,0elsewhere.
ORBIT208nonzeroqpatterns,remaining35fixedsymmetrypatterns vanish.
Projectedfunctiondimension (243+1-27-9)/4=52. Exactbasis_identity.json
reconstructsall52chargepolynomials onall1024physicalstates; span_rank.json
proves9sparseshapes+52completechargecolumns rank61 (all36newindependent).

Energyv10 signed_charge_telescope field bounded52canonicalexactnonzeroterms;
oldv9indicator/v8pair/v7quadfieldsremain. The new chart fixesold diagonal
profiles/compactcoefficients becausealltheirvariationslieinfull52span.
Twohoppingprofiles,twopenalties,ninesparsecoefficients,52patterncoefficients
optimized=66variablesincludingell. Allfourruns240iterations(hitlimit),under
hard250matrixevalcap;44newpatterncoefficientsnonzero. FreshCARaffinediff
<=3.56e-15. NOconvergence/controlledconditioningorperformanceclaim.

NewFAMILYv6 checksall52patternmomentszero plusallpreviousphysicalconstraints.
Cap66mixturesources (30old+36new). Actualsources59,58,59,58. ExactPERIODIC
lower→familyceiling gaps:
 W0 3.905668618334712e-5; W+.1 4.256669255112121e-5;
 W-.1 5.164377015198710e-5; W1 8.273009059499841e-5.
These areLOOSE comparedpriorfamilies; doNOTclaimcurrentloweristightlyoptimized
orfullsignedfamilyexhausted. Fixedsources/ratio/ceilings/sparse9;all52charge
coefficients/reflectedprofilescovered;notphysicalenergyuppersorotherfamilies.

Initialnear-eigenvectorLPpoolsfailedinfeasibleforall4; initial_dual_diagnostic
preserved. Addedall4096physicaldeterminantanchors as boundedfallback, producing
4548candidates. Physicalvacuumanchor guaranteesfeasibility. Exactrectangular
basis solverhandlesfewerindependentcolumnswithoutdroppinganyexactcondition.
Actual59/58basisvectors fromnonsymmetricphysicalstatesareallowed. Independent
productionreplayverifiesallmoments,energies,fidelities,trace/PSD/sourcebounds.
Newdriver signed_charge_family_numeric.py --determinant-anchors, cap5296
candidates; max66exactbasisrows. Scalarproposalworkisnotanacceptancebypass.

CONICFAILURE: oneW0CVXPY/CLARABELactive-conecheckattempted (settings2rounds,
40solverseconds/60iters/max16blocks). Lastobservedlive285swall; nativeprocess
exited137withoutanyproposal/logoutput. CauseNOTestablished. Agentdidnotkill;
processabsentwhenwallguardchecked. Noexactresult/convergence/improvement.
W_zero/conic_attempt_control.json recordsfailure; conic.logempty. Do not
restartfromobservationtimeout; thisprocessisconfirmedterminal. Solverneeds
boundeddiagnosis/replacementifpursued, notclaimedasanoptimizationproof.

EXACTCLASSICALEXTENSION: charge_markov_extension.json foreachcaseconstructs
PH+reflectionAVERAGED six-signed-chargeprobabilities. Rawprefix/suffixL1
mismatches .112145,.095655,.114194,.109537, soaveragingisexplicitlynecessary.
Afteraveragingall243prefix/suffixequalitieshold. Exact243-stateorder5Markov
kernelhasstochasticrows,shiftcompatibility,piP=pi,andrecovers729probabilities.
Positive5chargestates205,203,207,205;positive6patterns443,439,435,435.
Onlyclassicalaveragedchargesextend;fixedtotalN,spins,coherencesnotcertified.

NEXT CONCRETE OFF-DIAGONAL OBSTRUCTION (all4 exact coherent_overlap.json):
 B(i,j)=sum_spin(c_i^dagger*c_j+c_j^dagger*c_i).
 Five-siteY=B(0,3)-B(1,4); six-siteT=B(0,3)-2B(1,4)+B(2,5).
T12CARwords, normEXACT8 (trianglebound andinteger64component+8eigenvector).
Hermiticity,PH/reflection/spinflipinvariance checkedonall4096Fockstates;
periodictranslatedTsumzero. Moments W0 .00214708198962347,
W+.1 .0025722103258802843,W-.1 .0026907347054167144,
W1 .019802676676793113. ViolationpersistsunderPH/reflectionquantumaveraging.
Thusaveragedsignedchargesadmitclassicalstationaryextensionbutaveragedquantum
marginaldoesNOTadmitstationaryextensionmatchingthat6sitemarginal. Notaclaim
aboutinhomogeneousextensionofasinglemarginal/generalrepresentability.

Nextboundedattack: add compactmean-zero range-three HOPPING telescope
[gamma,-2gamma,gamma], withoutaddingphysicalrange3hoppingtotarget. Physical
symmetryblockcompatibilityalreadychecked, butproductionenergyintegration,
focusedtests,exactloweracceptanceandnewmatchedfamilymomentremainUNDONE.
Currentv6fullchargefamilyceilingswouldNOTcovernewoffdiagonalterm. One new
independenthoppingmomentcouldallow+1source. Otherconditionedhopping/spin
constraintsremainbeyondthissimpleoperator. Thecoarsecurrentfamilynumerical
limit isalsoanopenoptimizationtask; avoidblindlyclaimingscalarspaceexhausted.

Artifacts root results/marginal_graded_hubbard8/full_signed_charge/:
combined_summary.json, span_rank.json,basis_identity.json. EachW_zero,
W_plus_1_10,W_minus_1_10,W_plus_1 containsv10energycert,hints,range_two_replay,
v6familycert+range_two_family_limit_replay,charge_markov_extension,
coherent_overlap. No polishedsubdirs used; W_zero/conic hasNOproposal.
All120source/inputhashespercaseplus57basis/rankhasheschecked,total537.
Earlierchangedproofsources savedpre_signed_sources/manifest.json; pre-anchor
signed_charge_family_numeric_before_anchors.py also retained.

Code: newproductionhelper+energyv10+familyv6 withfull_signed_chargebool
requiringalloldermodes; adds36sourcesonlynewmode; histogramaccumulationchecks
all52patternmoments efficiently. Newnumericaldriverssigned_charge_numeric.py,
signed_charge_family_numeric.py, failedbounded signed_charge_conic.py.
Exactread-onlyprobes signed_charge_span_rank.py,signed_charge_basis_identity.py,
charge_markov_extension.py,coherent_overlap_probe.py. Existingmatchedfamily
replay nowrecognizesnewkind/refusesnewfieldwitholderfamilymode.

VALIDATION:596fulltests+96subtests pass454.40s. Focused20tests13.22s,
12family/correctiontests10.94s. Full_logfull_signed_charge/full_validation.log.
Centralresults/marginal_final_validation.json→full_signed_charge updated.
Alljobsconfirmedterminal. Generalrepresentability,arbitrarymolecular/longrange/
higherDtransfer,requestedaccuracyscalabilityremainunproved. Goalstaysactive.

## Latest accepted checkpoint — full charge indicators and classical extension

Read research/marginal_full_charge_indicators.md. Supersedes older checkpoints.
Goal ACTIVE. Previous turn made verified progress. All jobs complete; no GPU ran.
Never contact archived old task. No GPU spend this turn.

BEST final artifacts are in full_charge_indicators/<target>/polished/.
Parent target dirs contain earlier accepted energy plus FAILED initial dual
sampling diagnostics; use polished receipts for current bounds/family claims.

Exact million-site OPEN U4,t1,V1/2 intervals/site:
 W0 [-.6432237895841160,-.6106763470511881], exact lower
 -160805947396029/250000000000000; gain4.059518217284e-5 fromprecedingpairbound.
 W+.1 [-.6438715054058292,-.6114511605830021], exact lower
 -16096787635145729/25000000000000000; gain4.069229552820e-5.
 W-.1 [-.6428444456974698,-.6099015335193740], exact lower
 -3214222228487349/5000000000000000; gain4.005864816572e-5.
 W1 [-.6617035788877922,-.6184244823693281], exact lower
 -3308517894438961/5000000000000000; gain3.626556315036e-5.
All first+polished energy receipts accepted. All uppers freshlyrecomputed,
includingexact24siteenclosure, andunchanged. OriginalW0Hamiltonianunchanged.

New helper marginal_charge_indicator_telescope.py defines ALL_LABELS=12
canonical reflection-odd monomials of5binary p_i=q_i^2. LABELS contains6new
triples/quadruples:0,1,2;0,1,3;0,1,4;0,2,3;0,1,2,3;0,1,2,4.
EachY=productp_S-productp_reflectedS, localT=Y_left-Y_right.
Energyv9 uses higher_charge_indicator_telescope exactboundednonzero dict.
Older v8 squarepair andv7quadratic fields remain required. Sixnewindicators
complete12D indicator subspace (2singles+4pairs+4triples+2quadruples).
Exactrank12on32binarypatterns andrank25 foractual9sparse+6quadratic+4pairs+
6newdirections onall4096Fockdeterminants accepted in full_indicator_rank.json.
Sparse64entries unchanged, indices[1,2,4,5,6,7,31,11,20]. MaxPSD200,94local
blocks,total4096,151totalPSDchecks unchanged. No extraopening penalty.

New FAMILYv5 checks all12indicator momentszero andallpreviousconstraints;
30sourcecap vsold24, conditional full_charge_indicators boolrequirespairmode.
Allactualmixtures30states. Exact PERIODIC lower→fullindicatorfamilyceilings:
 W0 1.7942465185608435e-6; W+.1 1.8776663979155050e-6;
 W-.1 1.6584185985964500e-6; W1 7.714695899635311e-7.
Notphysicalenergyuppers. Fixedsources/ratio/ceilings/sparse9; unrestricted
indicator+quadratic coefficients/reflectedmean-correctprofiles arecovered.
Finiteceilingtolerance~2e-6, coarserthanprior1e-7. Exactoptimum/convergence
NOTproved. Othercorrectionspaces/sources/supports canescapetheseceilings.

NUMERICAL FAILURE/RECOVERY: first4optimizationshit120iterations, eachimproved.
Defaultdualradius.002 returned3infeasible+1unknown on336sampled states;
notmathematicalinfeasibility. initial_dual_diagnostic.json preserved.
Polishpass120iters forW0,+.1,-.1 (stillnotconverged),94forW1(converged).
Perpasshard250matrixevalcap unchanged; actualtotal240/target,214forW1.
Widerdualradius.02 produced30-vector exactrationalcandidates for all4.
Refinement.01for3old/.005forW1tested; selected .02old,.005W1 (lowercap).
Bothcandidateproposal/diagnostic files retained underpolished/. Fullindependent
physicalfamilyreplayacceptsselectedcaps. No furtherpolishing done.

NEW MATHEMATICAL RECEIPT: indicator_markov_extension.json ineachpolisheddir
proves exact stationary classical order5Markovextension ofbinaryp projection.
All32 prefix/suffixequalities, stochasticrows, shiftcompatibility, piP=pi,
andall64reconstructed6bitprobabilitieschecked. All32stationarystates and64
patternspositive. Threeindependentknownmodelchecks:iidbinary,alternating,
refusalofinconsistentdistribution. This extendsONLYbinaryindicatorlaw,
notqsigns,spins,fermioniccoherencesororiginalquantumdensitymatrix.

NEXT OBSTRUCTION: exact52PH-evenreflectionoddcharge-polynomialprobe finds
16zeromoments (12indicators+4otherquadratics),36nonzero oneachdual.
W0,+.1,-.1 strongestnormalized mixedsigned-charge separator:
 Y=q0*q4*(p2*p3-p1*p2); localnorm1,Ysupport64.
Moments .00014425139237642519,.00012921625861412112,.0001599940009566601.
W1 strongest Y=q0*q4*(p3-p1),localnorm2,Ysupport128,moment-.0003867305860386484.
Thusclassicalbinaryextensionexists butnecessarysignedchargeconsistencystill
fails: thesequantummixtures have no stationaryextension matchingtheir6site
marginal. Do NOTclaimnonextendibilityofasinglemarginalintoinhomogeneousstates.
No energyoptimization ofnewmixedsignedcharge operators yet.
Possible nextboundedattack: compact mixedcharge/indicator monomials; support64
foroneYdoesNOTmeanunionwithcurrent64sparseentryspanfitscap. Preserveexact
polynomial/PH/reflection/coefficient/PSD/refusalchecks. Fullchargepolynomial
space52 islargerthanindicator12; noncommuting/spin constraintsremainbeyondit.

Artifacts root results/marginal_graded_hubbard8/full_charge_indicators/.
combined_summary.json points toeach<case>/polished/. Thosecontainv9energycert,
hints,range_two_replay.json,v5familycert+range_two_family_limit_replay.json,
charge_polynomial_overlap.json,indicator_markov_extension.json.
Rootfull_indicator_rank.json. Source90hashentriesperfinaltarget+30rankhashes
checked; includingparentfirstenergyreceipts total566hashentriesverified.
Earlierchanged sourceversions savedpre_indicator_sources/; pre-widedual
samplingcopy diagonal_family_limit_numeric_before_wide_sampling.py also saved.

Productionenergyv9/helper, familyv5/sharedfull_charge_indicatorsbool.
Discoveryjoint_profile_numeric.py build addsindicator_terms, CLI
--full-charge-indicators requires--charge-square-pairs andalloldermodes.
Numericdual newfullindicatorflag/30rows and --perturbation-scale in(0,.2],
default.002. Matchedfamilydriverrefusesnewenergyfieldwitholderfamilymode.
Exactread-onlydiscoveryfull_indicator_rank.py andindicator_markov_extension.py.

VALIDATION:584fulltests+96subtests passed501.18s. Focused22energytests9.84s,
9familytests10.89s. Full logfull_charge_indicators/full_validation.log.
Centralresults/marginal_final_validation.json→full_charge_indicators updated.
All jobs finished. Generalrepresentability,arbitrarymolecular/long-range/higherD
transfer, andrequestedaccuracyscalabilityremainunproved. Goal staysactive.

## Latest accepted checkpoint — charge-square pairs and W=1 transfer

Read research/marginal_charge_square_pairs.md. Supersedes older checkpoints.
Goal remains ACTIVE. Previous turn made verified progress. No GPU ran/spend.
Never contact archived old task.

Best exact million-site OPEN intervals/site U4,t1,V1/2:
 W0 [-.6432643847662889,-.6106763470511881], exact lower
 -16081609619157221/25000000000000000; gain9.027242567648e-5.
 W+.1 [-.6439121977013573,-.6114511605830021], exact lower
 -8048902471266967/12500000000000000; gain9.902535166736e-5.
 W-.1 [-.6428845043456355,-.6099015335193740], exact lower
 -2009014076080111/3125000000000000; gain8.238359661880e-5.
 NEW W+1 [-.6617398444509426,-.6184244823693281], exact lower
 -4135874027818391/6250000000000000; width.043315362081614406.
All4 exact lower+upper accepted. Physicaluppers freshly computed with24-site
exact contraction inside160bit enclosure beforemillion recurrence. Threeold
uppers unchanged. W1 uses sameprojectors/sparse span/compactoperators/trialstate
recipe; lower scalarretuned. Comparedwith historicalnormperturbative width
2.0328451372 closes97.8692%; this is a looseboundcomparison, noexactsolution.
Accuracydegrades relativeW0 width.03258803771510083; notautomatictransferprecision.

New v8 energy field charge_square_pair_telescope supportsfour canonicalkeys
0,1;0,2;0,3;1,2. Define p_i=q_i^2 (empty/double indicator).
 Y_ij=p_i*p_j-p_(4-j)*p_(4-i), localT=Y_left-Y_right.
Eightproducts total; eachYsupport384,256,384,256 andeachlocalnorm2.
All translated corrections exactlycancel, originalW0Hamiltonianunchanged.
New helper experiments/marginal_charge_square_pairs.py reuses quadratic key/
coefficient validation, rejects single-sitepairs/malformed/nozero terms.
Exactrank10 ofsixquadratic+fournewpairs onall4096determinants verified.
SparseY still9shapes/64entries indices[1,2,4,5,6,7,31,11,20].
MaxlocalPSD200,94localblocks,total4096,151totalPSDchecks; no capincrease.

New v4 range2 FAMILY validatesall4squarepairmomentszero inadditiontosixquad
andNN/range2profiles. Mixturecap24 vsold20, eachactualdual24integerstates.
All4matched exactfamilycaps accepted; PERIODIClower→familyceiling gaps:
 W0 1.1348724240303176e-7; W+.1 8.162929920975288e-8;
 W-.1 9.313586930791479e-8; W1 8.777598536513510e-8.
These cover unrestrictedcoefficients forfullquadratic+all4squarepairs,
reflectedmean-correctNN/range2profiles, fixedprojectors/ratio/ceilings/sparse9.
Notphysicalenergyuppers or otherfamily/support limits. Olderfamilyversions
keepoldcaps; matchedreplayrejectsenergy squarefieldswithoutfullsquaremode.

NEXT: exact52-polynomialprobe nowfinds10zeromoments,42nonzero onall4duals.
Strongest normalizedseparator sameforall:
 Y=p2*p3*p4-p0*p1*p2, localnorm2,Ysupport192,two tripleproducts.
Moments W0 .0011070092603546358,W+.1 .0011247038548555521,
W-.1 .0010830935598907289,W1 .0010517597452266655.
No triple optimization/energyimprovement tested.
indicator_basis_summary.json selectsfull12-dimensional PH-evenreflectionodd
functions of5binaryp_i:2singles+4pairs+4triples+2quadruples.
Singlescoveredbyquadraticprofiles,pairsnowcovered; remaining6allnonzero.
Concrete nextboundedattack: testallremaining6indicator directions together,
then mixed signedcharge/noncommutingconstraintsstillremain. Couldgeneralize
compactindicatorformula beyondpairs usingboundedcanonical subsets, preserving
exactkey/degree/source/PSD/refusalchecks. Do notclaimgeneralrepresentability.

Artifacts root results/marginal_graded_hubbard8/charge_square_pairs/:
combined_summary.json;charge_square_pair_rank.json;indicator_basis_summary.json;
subdirs W_zero,W_plus_1_10,W_minus_1_10,W_plus_1 eachenergycert,hints,
range_two_replay.json,range_two_family_limit_certificate.json,
range_two_family_limit_replay.json,charge_polynomial_overlap.json.
All85source/inputhashes pertarget checked;rankreceipt28hasheschecked.
Historicalchanged sources savedpre_square_pair_sources/manifest.json.
Production: energyv8; familyv4 shared_replay charge_square_pairs boolrequires
full_quadraticmode, adds4momentsand4sources onlyinnewmode.
Discovery joint_profile_numeric.py build accepts square_terms and new
--charge-square-pairs (requires--quadratic-charge andfree-range2/W).
Numericdual sameflag,24rowLP/distinctuntrustedproposaltype; allproposalbasis
weights solvedexactly, thenindependentphysicalreplay. Numericsearches
115,115,114,119evaluations under250cap, freshCARaffineerror<=9.06e-14.
Production source code stable; no edits during proofreceipt hashing.

Focused22energy/correctiontests6.34s and9familytests11.43s passed.
Full regression passed570tests+96subtests in416.07s. Log:
results/marginal_graded_hubbard8/charge_square_pairs/full_validation.log.
Allresearchjobscomplete; no localcompute running. Centralvalidation updated
to charge_square_pairs. Final372source/inputhashentries checked.
General molecular/arbitrarylong-range/higherD transfer,representability,
andrequestedaccuracyscalabilityremainunproved.

## Latest accepted checkpoint — compact quadratic charge telescope

Read research/marginal_quadratic_charge_telescope.md. Supersedes older checkpoints.
Goal remains ACTIVE; this turn makes concrete progress, does not solve broad goal.
Never contact archived old task. No GPU ran/spend this turn.

Best exact million-site OPEN intervals/site, U4,t1,V1/2:
* W0: [-.6433546571919653,-.6106763470511881], lower exact
  -16083866429799133/25000000000000000; gain7.534081772296e-5.
* W+.1: [-.6440112230530247,-.6114511605830021], lower exact
  -8050140288162809/12500000000000000; gain5.589655587208e-5.
* W-.1: [-.6429668879422543,-.6099015335193740], lower exact
  -8037086099278179/12500000000000000; gain9.607866993368e-5.
All lower+physical-upper replays accepted. Each upper freshly recomputed,
including exact24-site enclosure, and unchanged. W0 physical Hamiltonian unchanged.

v7 energy field quadratic_charge_telescope {'0,3': gamma} adds
 Y=q0q3-q1q4, T=Y_left-Y_right=q0q3-2q1q4+q2q5.
Gammas W0 -492/3125, W+ -66823/500000, W- -180299/1000000.
New helper experiments/marginal_quadratic_charge_telescope.py supports all6
canonical reflection-odd quadratic pairs (0,0),(0,1),(0,2),(0,3),(1,1),(1,2).
Y has416 entries but2 products. SparseY stays9shapes/64 entries with indices
[1,2,4,5,6,7,31,11,20]. MaxPSD200,94local blocks,total4096;151allPSDchecks.
No opening penalty for telescope; same physical2t+|V|+2|W| deduction/N.

Exact v3 range-two FAMILY caps now check all6 quadratic momentszero.
Five are implied by profiles; exactly one extra independent row (20 total).
All20-vector rational mixtures accepted. Remaining PERIODIC lower→familycap:
 W0 5.406925868495681e-8; W+ 6.144059153559798e-8;
 W- 1.233973613464014e-7.
Caps cover all6 quadratic coefficients, all reflected mean-correct NN/range2
profiles, fixed sources/ratio/ceilings and9sparse span. Not physicalenergyuppers.
Do not apply v1/v2 familycaps to new energyquadfield; matched replay refuses.

Exact all4096-CAR identity replay quadratic_profile_identities.json proves
quadratic telescope rank6, existingprofile spanrank5: labels 00,01,02,03,11,12
have identities 2(A-B),D-E,R2,R3,2B,E. A,B onsite, D,E NN density derivatives.
Receipt32sourcehashes checked. Every energy/family/separator case83hashentries
checked. Historical sources saved pre_quadratic_sources/manifest.json.

NEXT EXACT SEPARATOR: full52-element PH-even reflection-odd five-site charge
polynomial basis with exponents0,1,2. All6 quadratic momentszero; allother46
nonzero on each accepteddual. Strongest normalized moment:
 Y=q3^2*q4^2-q0^2*q1^2; T=Y_left-Y_right; norm2, Ysupport384,two products.
Moments W0 -.00355248501238245, W+ -.003784044864915257,
W- -.00333752478336591. q^2 indicates empty/doublyoccupiedsite.
Exact probe independently checks expansion/reflection/PH overFockspace.
No quartic optimization or energy improvement tested yet. Next boundedattack
could compactly add this degree4 charge consistency direction; preservecaps
and representability refusal paths. Other noncommuting/spin constraints remain.

Production: energy v7 extends v6 with bounded compact term dict, allFockchecks.
Shared family _replay(...,full_quadratic_charge=False) only permits newmode
with free_range_two_profile; checks all6moments, one extra source20 vs19.
New family v3 invokesmode; older versions keep previous caps.
Numeric joint_profile_numeric.py --quadratic-charge requires explicitW and
--free-range-two-profile, addsD3 scalar; build now accepts quadratic_terms.
Numeric diagonal_family_limit_numeric.py sameflag adds20thLProw and distinct
untrusted joint_diagonal_quadratic_charge_family_proposal_v1.
Existing range_two_density_replay.py permits v7. Existing matchedfamilydriver
recognizes newproposal→v3 and requires fullquadcap forenergyquadfield.

Artifacts root results/marginal_graded_hubbard8/quadratic_charge_telescope/.
combined_summary.json compares all3. W_zero/,W_plus_1_10/,W_minus_1_10/ each
contains profile_joint_r1_2_certificate.json, congruence_witnesses.json,
range_two_replay.json, range_two_family_limit_certificate.json,
range_two_family_limit_replay.json, charge_polynomial_overlap.json.
Discovery scripts charge_polynomial_overlap_probe.py and
quadratic_profile_identities.py are exact read-only probes, no acceptance bypass.

Focused22correction/range2 tests passed9.95s;22family/correction tests18.37s.
Full regression passed556tests+96subtests in386.06s. Log:
results/marginal_graded_hubbard8/quadratic_charge_telescope/full_validation.log.
No local research jobs remain. Central validation updated→quadratic_charge_telescope.
General molecular/arbitrary long-range/higher-D transfer, representability,
and requested-accuracy scalability remain unproved.

## Latest accepted checkpoint — free range-two profiles

Read research/marginal_free_range_two_profile.md. This supersedes older checkpoints.
Goal remains active. Previous turn made progress. No jobs or GPU running.
Never contact the archived old task. No GPU spend this turn.

CURRENT best million-site open intervals/site (U4,t1,V1/2):
* W=0: [-.6434299980096883,-.6106763470511881], lower exact
  -16085749950242207/25000000000000000. Original nearest Hamiltonian unchanged.
  Gain9.548622846432e-5 over preceding lower; closes .290681% preceding width.
* W=+1/10: [-.6440671196088968,-.6114511605830021], lower exact
  -805083899511121/1250000000000000. Gain .00013695253351456.
* W=-1/10: [-.6430629666121880,-.6099015335193740], lower exact
  -160765741653047/250000000000000. Gain5.998238944816e-5.
All physical uppers recomputed with exact24-site enclosure and million-site
transfer; same as previous uppers. All actual lower+upper replays accepted.

Added freedom: local profile [5W/4+delta,5W/4-delta,5W/4-delta,5W/4+delta].
At W0 delta=-39877/250000: the entire added profile sums to zero globally.
It is the compact telescope Y_left-Y_right, Y=q0*q2-q2*q4 (320 diagonal
entries but2 charge products). The v6 lower format already supported the
signed profile; lower production and PSD caps did not need changes.
Sparse-Y remains64 entries, max local PSD200,94local blocks/151total checks.

New v2 range-two FAMILY verifier allows free profiles and requires exact
zero expectation of [1,-1,-1,1]. With9 shapes,19scalar constraints and<=19
physical mixture sources; old fixed-profile version remains<=18. All19-state
mixtures accepted exactly, including trace/fidelity/six nearest gradients,
nine sparse moments and the additional free-profile moment.
Remaining periodic gaps to fixed-source/span FREE-PROFILE family cap:
 W0 9.794610292439077e-8; W+ 7.815154859454668e-8; W- 1.1177187240864765e-7.
These are caps on LOWER certificates, not physical energy uppers or limits
for other supports/sources. Existing v1 fixed-profile caps remain valid only
for their smaller families.

NEXT CONCRETE COMPACT SEPARATOR: exact probe of all6 reflection-odd quadratic
charge polynomials on5sites finds5zero moments and1nonzero:
 Y=q0*q3-q1*q4; T=Y_left-Y_right=q0*q3-2*q1*q4+q2*q5.
Moments W0 -.004784401291454263, W+ -.0041811892879441166,
W- -.00533049067134383. Y has416 nonzero diagonal entries but only2 products.
No optimized bound using this direction yet. Could explore a compact CAR/
charge-polynomial telescope representation or a mean-zero range-three local
profile; do not blindly raise sparse support or PSD caps. Scalar tuning in
current family is now capped. Broader transfer and general claims still open.

Validation:542 full tests+96subtests pass392.73s;19focused pass13.85s.
All energy/family/quadratic-separator source hashes checked,81 entries/target.
Full current suite includes earlier supplemental tests. Old executed sources
are saved under free_range_two_profile/pre_free_profile_sources/.

Artifacts: results/marginal_graded_hubbard8/free_range_two_profile/combined_summary.json.
Subdirs W_zero/, W_plus_1_10/, W_minus_1_10/: profile_joint_r1_2_certificate.json,
congruence_witnesses.json, range_two_replay.json,
range_two_family_limit_certificate.json, range_two_family_limit_replay.json,
quadratic_charge_overlap.json, combined_summary.json.
Central status results/marginal_final_validation.json -> free_range_two_profile.

Production changes this turn: marginal_diagonal_family_limit.py public replay
keeps old bounds; private shared _replay can check one physically specified
extra free-profile moment and allows one extra mixture source only in that
mode. marginal_range_two_family_limit.py v2 calls that mode; v1 remainsfixed.
Discovery joint_profile_numeric.py --free-range-two-profile addsdelta and
records it; build now accepts optional exact range_profile. 250eval cap and
fresh physical reconstruction retained. diagonal_family_limit_numeric.py
--range-two --free-range-two-profile adds the19th moment/LP row and distinct
untrusted proposal kind. range_two_family_limit_replay.py recognizes it and
uses accepting v2. New quadratic_charge_overlap_probe.py checks compactYs.
No local runtime jobs left. General molecular/arbitrary long-range/higher-D
transfer, representability and requested-accuracy scalability remain unproved.


## Latest accepted checkpoint — range-two density transfer

Read research/marginal_range_two_transfer.md. This supersedes older checkpoints.
Goal remains active. Previous turn made progress. No local job or GPU running.
Do not contact the archived old task in either direction. No GPU spend this turn.

Two NEW Hamiltonians: open half-filled U4,t1,V1/2 chains with W*q_i*q_(i+2).
Same original projector sources, nine selected shapes, physical block state and
boundary filter recipe; scalar lower profiles/penalties retuned.
* W=+1/10 million-site interval/site:
  [-.6442040721424114,-.6114511605830021], width .03275291155940926.
* W=-1/10 interval/site:
  [-.6431229490016361,-.6099015335193740], width .03322141548226207.
Both have fresh exact lower AND physical upper replay, including exact24-site
contraction enclosed before the million-site evaluation. ~22.2s each.
They close85.93%/85.73% of the conservative norm-perturbed old interval width.
This is range-two density transfer within1D, not arbitrary long-range/molecular.

Fixed UNIFORM range-two profile family caps also accepted: 18-state PSD dual,
all six nearest-profile and nine diagonal moments zero, fidelity/trace exact.
Remaining periodic family gaps9.446437040882777e-8 (W+) and5.122156012051345e-8 (W-).
The range-two profile is FIXED to [5W/4]*4 in these caps. They are upper limits
on lower certificates, NOT physical energy uppers or free-profile-family caps.

NEXT CONCRETE LEAD: free range-two profile [lambda,5W/2-lambda,5W/2-lambda,lambda].
The accepted fixed-profile mixtures have exact derivative expectations
-.007055426323878254 and-.004769296896628078 for direction [1,-1,-1,1].
Exact4096-state check identifies it as Y_left-Y_right with Y=q0*q2-q2*q4.
This Y has320 diagonal entries but only2 charge-product terms. Production v6
already accepts the signed mean-correct4-entry range-two profile, so a new
numeric free-profile parameter can test this compact constraint without raising
any sparse-Y or PSD cap. It can also be explored at W=0 on the original target.
Free-profile numerical optimization and its19-constraint family dual are NOT
implemented/accepted yet. Do not conflate the current fixed-profile cap with it.

Production:
* marginal_projector_extendibility.py accepts explicit v6 with local kind
  local_hubbard_range2_block_v1, target W and optional range_two_density_profile.
  Exact reflected profile sum5W at L6; local all-Fock diagonal addition.
  Opening norm2t+abs(V)+2abs(W); old versions refuse range-two fields.
* marginal_range_two_density.py supplies local_profile and diagonal_value.
* marginal_range_two_transfer.py composes existing physical nearest compiler
  with six-site dressed range-two patches. Existing48D recurrence reused;
 4-site partial traces on each block; exact physical expanded tests4/8/12sites.
* marginal_range_two_family_limit.py wraps existing nearest diagonal-family
  constraint replay, adding exact expectation of a FIXED range-two profile.
  Old marginal_diagonal_family_limit.py now refuses W/profile fields; wrapper
  strips them before reusing nearest constraints and checks them separately.

Validation:534 tests+96 subtests pass368.15s. Subsequent19 family tests+12
subtests pass14.01s, covering5 new cases and final nearest-only field guard.
All actual energy/family/separator source hashes checked,81 entries per target.
Earlier executed source versions are saved in range_two_density/pre_range_two_sources/.
Central status results/marginal_final_validation.json -> range_two_density.

Artifacts: results/marginal_graded_hubbard8/range_two_density/combined_summary.json.
Subdirectories W_plus_1_10/ and W_minus_1_10/ have v6 certificate, optional
congruence_witnesses.json, range_two_replay.json, range_two_family_limit_certificate.json,
range_two_family_limit_replay.json and free_profile_separator.json.
Drivers in discovery/: joint_profile_numeric.py accepts --target-w (currently
uniform local range-two profile). diagonal_family_limit_numeric.py requires
--range-two for these targets and produces a distinct untrusted range2 kind.
range_two_density_replay.py, range_two_family_limit_replay.py and
range_two_profile_separator.py perform separate standard-library exact checks.
Old nearest-only family replay driver explicitly refuses range-two targets.

Original W=0 best remains -.6435254842381526 per site, historical upper
-.6106763470511881. Previous negative-V transfer remains accepted. Generic
molecular/arbitrary long-range/higher-D transfer, global representability and
requested-accuracy scalability remain unproved. Keep broad goal active.


## Latest accepted checkpoint — support replacement and transfer

This supersedes historical checkpoints below. Read research/marginal_shape_selection.md.
Previous turn was progress. Goal remains active; no local jobs or GPU running.
Do not contact the archived old task in either direction.

* Current best million-site U4,t1,V1/2 lower/site:
  -3217627421190763/5000000000000000 = -.6435254842381526.
  Historical upper -.6106763470511881; width .0328491371869645.
  This turn gained .00013407268524356 per site and closed .40649% previous gap.
* Best indices [1,2,4,5,6,7,31,11,20], nine shapes, exactly64 nonzero entries.
  No production modules or matrix caps changed. Five bounded CPU trials;
  two best stages received full exact energy and family-limit acceptance.
* Corrected previous bookkeeping: original nine-shape extra31 certificate
  had60 entries, not64. Its actual certificate and exact receipts were always
  correct. Shape31 has4 entries. New selected support really uses64.
* New fixed-family limit -.6435228599405247 PERIODIC/site, gap
  1.242976278974152e-7 above accepted periodic lower -.6435229842381526.
  Exact18-vector PSD mixture: all six profile and nine diagonal moments zero,
  trace1, fidelity inequalities checked. Not a physical energy upper.
* New mixture still violates93 of120 directions; top indices0,12,8 with
  moments approx .00192879,-.00185326,-.00156874. Exact separators are saved.
  Scalar tuning inside this family is now capped. Further gains need other
  supports, sources/ratio, or a different bounded correction representation.
* Same nine shapes/sources transfer to U4,t1,V=-1/2 after scalar retuning:
  [-.5510143532722551,-.5240746138143995] per site. Improves previous accepted
  eight-shape lower by6.481601897812e-5, closes .24002% of its interval.
  Fresh physical upper replay includes exact24-site contraction enclosure.
* Exact matched energy15.478s; exact limit1.961s; exact transfer18.050s.
  Optional integer congruence witnesses, standard-library-only acceptance.
  36 focused tests+8 subtests pass6.70s. Earlier504-test regression+14-test
  supplement still apply: production modules unchanged this turn.

All new artifacts are under
results/marginal_graded_hubbard8/joint_projector/signed_density/shape_selection/.
selection_summary.json lists five trials and their acceptance status.
Best directory swap_3_for_20/ contains profile_joint_r1_2_certificate.json,
profile_joint_r1_2_certificate_accelerated_replay.json, congruence_witnesses.json,
diagonal_family_limit_certificate.json, diagonal_family_limit_replay.json,
remaining_diagonal_overlap_moments.json and combined_summary.json.
transfer_V-1_2/ underneath has matched_transfer_replay.json and combined_summary.json.
Intermediate exact stage is swap_0_for_11/. Other3 trials are numerical only.
Current central validation: results/marginal_final_validation.json -> shape_selection.

Discovery joint_profile_numeric.py now accepts --shape-indices, --seed-certificate,
--output-dir. It enforces distinct bounded indices and a<=64-entry union.
The selected indices and seed hash are recorded. joint_energy_replay.py and
joint_transfer_replay.py accept optional --psd-witnesses JSON; all exact gates
remain in production replay. Energy driver writes *_accelerated_replay.json
for optional witnesses so an original elimination receipt is never overwritten.
Three diagonal-family drivers now accept --directory; numeric requires9 shapes
and obtains them from profile_proposals.json. Replay matches source/ratio/
ceilings/span to accepted energy. Old executed driver snapshots and hashes
are saved under shape_selection/pre_generalized_dual_sources/ and pre_*.py.

Generic molecular, long-range and higher-dimensional transfer, global
representability and requested-accuracy scalability remain open. Only1D
nearest-neighbor parameter transfer has been shown. Avoid conflating finite
families, fixed-size proof speed, and general quantum chemistry. No new GPU
spend this turn; previous instance remains terminated, historical ~$0.72 estimate.


## Latest accepted checkpoint — nine-shape limit and fast exact replay

This block supersedes all historical checkpoints below. Goal remains active.
Read research/marginal_nine_shape_limit.md and research/marginal_exact_congruence.md.
No GPU is running. Do not message the archived old task in either direction.
No local proof/search/test job remains running at this checkpoint.

* Accepted million-site U4,t1,V1/2 lower/site:
  -2011436115385613/3125000000000000 = -.6436595569233962.
  Historical physical upper -.6106763470511881; width .03298320987220806.
  Nine diagonal shapes,60 nonzero entries, within the existing64-entry cap.
* Full original exact replay accepted in401.725s. Optional integer congruence
  replay accepted all151 blocks in15.498s, same bound and all ranks/nullities.
  73 unique factors used79 times;72 fallback blocks. Auxiliary data4.08MB,
  energy certificate18,879 bytes, proposal7.960s. Observational runtime ratio
  25.92, not repeated controlled benchmark or scaling theorem. Stdlib acceptance.
* New exact18-state PSD mixture caps the FIXED nine-shape family at
  -.6436569911639332 per periodic site, leaving only6.575946298619209e-8
  above accepted PERIODIC lower -.6436570569233962. This is an upper limit
  on lower certificates, NOT a physical ground-energy upper. Every one of
  six profile gradients and nine diagonal moments is exactly zero; trace1
  and both fidelity inequalities exact. Sources/ratio/ceilings/span matched.
* This dual still violates93 of120 other tested consistency directions.
  Largest is basis index11, exact moment approx .0039845523827132835 with
  eight entries. Current64-entry cap prevents straightforward addition.
  Next useful work: bounded shape replacement/compression or different source/
  ratio families, broader transfer; do not continue scalar tuning inside the
  now-capped family. No stronger bound from these remaining separators yet.
* Regression504 tests+96 subtests passed in349.66s. New diagonal-family
  supplemental14 tests (11 new,3 existing)+12 subtests passed in7.72s.
  Existing project venv stalled on SciPy library reads; those processes were
  terminated and jobs rerun using system python with OPENBLAS_NUM_THREADS=1.
  System python is /opt/homebrew/Caskroom/miniconda/base/bin/python (3.12.2).
  No test failures or lingering stalled jobs remain.

All nine-shape receipts are under
results/marginal_graded_hubbard8/joint_projector/signed_density/symmetric_diagonals_9/extra_31/:
profile_joint_r1_2_certificate_replay.json, congruence_witnesses_replay.json,
congruence_summary.json, diagonal_family_limit_certificate.json,
diagonal_family_limit_replay.json, remaining_diagonal_overlap_moments.json.
Modified pre-congruence sources retained in joint_projector/pre_congruence_sources/.
Current central status: results/marginal_final_validation.json.

New production modules: experiments/marginal_congruence_psd.py,
experiments/marginal_diagonal_family_limit.py. The latter reuses old mixture
and shape validators and allows9+k mixture sources for at most9 shapes;
old v1 family verifier's16-source cap is unchanged. Full-Fock PSD caps unchanged.
New optional psd_witnesses keyword on extended energy replay retains original
Bareiss fallback and all singular/exact-division refusals. New family-limit
module requires exact zero expectations for all supplied bounded-support shapes.

Discovery drivers: joint_congruence_replay.py (proposal vs stdlib replay),
diagonal_family_limit_numeric.py (bounded244-state/18-basis proposal),
diagonal_family_limit_replay.py (exact physical independent replay),
diagonal_family_overlap_probe.py (exact remaining separators), all in
results/marginal_graded_hubbard8/discovery/.

Earlier exact negative-V transfer still valid: same eight shapes/sources,
[-.5510791692912332,-.5240746138143995], closes37.53845% previous target gap.
Only1D nearest-neighbor parameter transfer is established. Generic molecular,
long-range, higher-dimensional transfer, representability, and requested-
accuracy scalability remain unproved. Do not mark broad goal complete.


## Latest checkpoint — diagonal corrections and completed GPU run

This supersedes the older checkpoint below. Read
`research/marginal_diagonal_overlap.md` and `research/marginal_gpu_diagonal_search.md`.
The goal remains active and full generality/scalability are unproved.

Accepted results:

* New exact U4,t1,V1/2 million-site lower/site:
  -201152653303269/312500000000000 = -.6436884905704608.
  Historical physical upper -.6106763470511881. Eight diagonal correction
  shapes,56 nonzero entries, cross the preceding family limit by
  .00030075000678818185 per PERIODIC site. Exact receipt:
  joint_projector/signed_density/symmetric_diagonals_8/profile_joint_r1_2_certificate_replay.json
  (relative to results/marginal_graded_hubbard8/).
* New v5 accepting path in experiments/marginal_projector_extendibility.py
  checks reflection-odd five-site diagonals with at most64 entries, adds
  exact T=Y_left-Y_right within all94 local blocks, max200. No matrix caps
  raised. Three new tests include every eight-site determinant's cyclic
  cancellation. All482 tests passed in485.605s;13 focused tests in10.437s.
* Transfer accepted with SAME projector sources and eight shapes at
  U4,t1,V=-1/2: interval/site [-.5510791692912332,-.5240746138143995].
  Fresh exact lower and physical upper replay;37.53845% narrower than this
  target's previous interval. Receipt:
  symmetric_diagonals_8/transfer_V-1_2/matched_transfer_replay.json.
  This is one-dimensional parameter transfer, not generic molecular transfer.
* Exact Gram congruence performance lead: one actual132x132 block accepted
  via nonsingular integer triangular R and strict diagonal dominance of
  R^T A R. Integer check .461s; independent python -S reconstruction and
  generic multiply3.893s total. See joint_projector/congruence_probe/.
  Not yet integrated into full energy replay; preserve singular PSD fallback
  and refusal tests if developing it. No speedup/scaling guarantee asserted.

GPU authorization and completed work:

The old task sent its FINAL message: user approved responsible Lambda
spending up to$2600 total and asked no further messages either way. That
old task is archived. DO NOT SEND IT PROMPTS, REPORTS OR ACKNOWLEDGMENTS.
A running A100 at129.146.98.55 was handed here. This task rebuilt current
physical tensors and actually evaluated785 LP-directed extra-diagonal
candidates across all94 sectors on A100 in2.8596s.25 CPU checks included
winners; max discrepancy7.501e-13. All downloaded hashes matched.
Results: results/lambda_runs/diagonal_science/.

The GPU instance has since been TERMINATED IN LAMBDA. The UI was checked
and showed No running instances. Observed lifetime estimate$0.72 at$1.99/h
includes the earlier benchmark period, not a provider invoice. No GPU is
left running. The existing SSH key remains usable for a future bounded
launch; do not assume the old IP is active. See migration/GPU_ACCESS.md and
results/lambda_runs/diagonal_science/lifecycle.json.

COMPLETED NINE-SHAPE EXACT JOB:

GPU screening selected additional basis shape31. CPU refinement of the
resulting nine-shape family proposed periodic density-.6436570569233961,
with64 nonzero diagonal entries and fresh CAR mismatch below9e-15.
Its standard-library exact energy replay ACCEPTED in401.724721458 seconds.
Million-site open lower -2011436115385613/3125000000000000
= -.6436595569233962, with historical upper -.6106763470511881.
Receipt: joint_projector/signed_density/symmetric_diagonals_9/extra_31/profile_joint_r1_2_certificate_replay.json.
All31 input/source hashes checked. Subsequent proof-module changes have
pre-change copies and manifest in joint_projector/pre_congruence_sources/.
The earlier PID59824/session85226 job is finished; do not restart it.

CURRENT DEVELOPMENT: optional exact congruence SPD witnesses are integrated
into the extended replay via keyword psd_witnesses. Original no-witness
Bareiss remains unchanged, including null pivots and exact divisions.
New experiments/marginal_congruence_psd.py accepts only bounded invertible
upper-triangular integer factors with exact strict positive diagonal
dominance, keyed to freshly reconstructed matrix SHA256. Unused hints refuse.
22 focused tests pass, including mixed singular-fallback and energy replay.
Full current regression and all-block accelerated replay still pending.
Discovery joint_congruence_replay.py proposes hints in an explicitly
NONACCEPTING mode with a subprocess-local matrix collection hook; its
replay mode is standard-library-only and never installs that hook.
Numerical witnesses generated for nine-shape certificate in7.960 seconds,
151 blocks collected, no Cholesky proposal failures. Exact full accelerated
replay is now launched; inspect congruence_exact_replay.log and actual
process before duplicating work. No GPU running.

Numerical script now supports --signed-density --symmetric-diagonals8
--extra-diagonal31 (spaces between option and integer in actual shell) and
--target-v=-1/2 for the separate transfer probe. GPU scripts live beside
results/lambda_runs/diagonal_science/science.npz and receipts. Main proof
module has optional congruence changes after the482-test validation. Exact rational CPU acceptance
remains mandatory.

## Earlier checkpoint (historical)

## Current accepted checkpoint — continued locally

Read `research/marginal_joint_projectors.md` and
`results/marginal_graded_hubbard8/joint_projector/signed_density/combined_summary.json`.
They supersede the older scientific status below and in CLOUD_HANDOFF.md.
The full objective remains active, not complete. No paid cloud resources or
archive upload occurred. All research compute processes have finished at this
checkpoint; do not restart from old running logs.

* New exact million-site open lower/site: -1006233292669113/1562500000000000
  = -0.6439893073082323. Historical physical upper remains -0.6106763470511881.
  Width0.03331296025704422; 6.1323% narrower than the preceding interval.
* Exact combined path `hubbard_projector_extension_v4` now exists in
  `experiments/marginal_projector_extendibility.py`. It combines residual
  half penalty alpha with joint beta*(Ph+rPc), retaining one rank-one update
  per spin sector. Joint r=.5 tight B=595736407/250000000 accepted in actual
  energy replay; all4096 states,94 local blocks,max200; joint1280cols,max132.
* Strongest certificate/replay:
  `joint_projector/signed_density/profile_joint_r1_2_certificate.json` and
  `_certificate_replay.json`, relative to results/marginal_graded_hubbard8/.
  Driver: discovery/joint_energy_replay.py (run under python -S).
* Numerical limit of THIS fixed-source,r=.5,fixed-ceiling profile family is
  now exactly bracketed within6.67309833381364e-8/site. Nine integer local
  vectors with exact positive mixture weights saturate both fidelity limits
  and cancel all six profile derivatives exactly. Limit upper is about
  -.643986740577249 PERIODIC; accepted periodic lower is -.6439868073082323.
  This caps attainable relaxation lower certificates, not physical energy.
  `experiments/marginal_joint_family_limit.py`, signed_density/family_limit_*
  and discovery/joint_family_limit_replay.py hold the accepted proof.
* Exact witness violates full five-site overlap consistency: diagonal entry
 102 of rho_left-rho_right is about-.008773379064595197. Reflection-odd
  Y=|102><102|-|612><612| gives telescoping expectation-.017546758129190394.
  Receipt: signed_density/overlap_consistency_separator.json. A single-Y
  SLSQP trial found zero coefficient/no gain; do not claim that extension
  works. Its proposal kind is deliberately unsupported by the verifier.
  Next scientific attacks: several simultaneous overlap constraints,
  changing r or source families, and transfer beyond the matched model.
* All479 tests passed in406.727seconds; complete log
  joint_family_limit_full_validation.log. Four exact receipts and123 source
  hashes checked. results/marginal_final_validation.json is updated.
* Correct root numerical script: discovery/joint_profile_numeric.py; supports
  signed density and exploratory one-diagonal telescoping. It uses analytic
  eigenvector gradients, step-normalized affine derivatives, and fresh CAR
  reconstruction. Agent full_profile_joint.py was WRONG (missing/.01),
  lacked its claimed fresh check, and reversed better/worse. See rejected
  diagnostics; never promote its outputs.
* Coarse independent joint overlap receipts at r=.5,1,2 accepted B2.384,3.04,
  4.926. Tighter r1 andr2 proposals remain unaccepted. Tight r.5 IS accepted
  by the current energy receipts. Do not confuse the two sets.

GPU work requested during continuation: bounded CuPy/CUDA benchmark prepared
in discovery/joint_gpu_benchmark.py; CPU reference validated, no GPU ran
(CuPy absent locally). Final CPU six-evaluation receipt: prep11.460854s,
scalar subset eigenvalues.105665s, batched.113940s, Gram.005810s. GPU harness
execution remains untested; no speedup claim. See research/marginal_joint_gpu_benchmark.md.
Existing Lambda SSH key/access was found by the old task; see
migration/GPU_ACCESS.md. Dashboard showed no running instance. User has about
$2600credit, which is not blanket spending authority. The old task is
coordinating a bounded instance launch scope; do not duplicate that flow.

General representability, requested-accuracy scalability, generic molecular
interactions and broader physical transfer remain unresolved.

## Original handoff (historical)

The user's latest instruction is to continue this goal and project in a forked or fresh LOCAL session. Cloud migration was canceled. A same-directory fork failed with a history projection ordinal error, so this checkpoint is for a new task in the existing Spectra folder.

The original task ID is 01a08e58-e30c-7800-8df6-ab713cbd7787. Its goal is currently paused. Do not resume competing writers in that old task. The new task should create an active goal with this full objective and continue working:

> Resolve the quantum-marginal certificate blockers: certify or refute the current relaxation's numerical limit, identify and test violated physical positivity constraints, seek compact stronger certificates, and evaluate transfer beyond the matched model. Track remaining unproved scalability and general representability claims explicitly; do not equate finite-model successes with solving general quantum chemistry.

Read `migration/CLOUD_HANDOFF.md` for the detailed mathematical checkpoint, evidence, current code, exact proofs and next attack. Its CLOUD filename and cloud setup instructions are historical: use this existing local directory and Python environment; no upload, VM, cloud purchase, or cloud Work task is wanted now. The scientific checkpoint remains current.

The user wants sustained original mathematical/computational work, multiple agents where useful, and concrete progress rather than literature-review assignments or repeated status summaries. The research representation under attack is the geometry of physical marginals. Be candid about what is proved. Preserve validation and refusal paths. Use only mattpocock skills unless the user changes that preference. Follow the guidance copied in `migration/agent-instructions/` when applicable.

Latest fully validated milestone: U4,t1,V1/2 million-site energy/site interval approximately [-.646165607326505,-.6106763470511881], narrowing the previous width by32.4402%. All469 tests passed at that milestone. Exact receipt: `results/marginal_graded_hubbard8/six_site_projector/refined_independent_replay.json`.

Subsequent work added `experiments/marginal_charged_projectors.py` and five focused tests in `tests/test_marginal_charged_projectors.py`, all passing in the final local run (0.296 seconds). No complete regression has run after that addition. The direct joint-bound test was corrected by root to retain the physical source norms, include both acceptance and refusal, and test amplitude-rescaling invariance. Do not revert it to the agent's earlier incorrect version.

Immediate next step: calculate and exactly replay a JOINT overlap bound for Phalf+r*Pcharged, then seek a stronger accepting energy certificate. The separate charged-penalty extension has a proved dominance obstruction: theta_half+theta_charged>1 and Phalf+Pcharged<=I. The full derivation and coupled alternative are in CLOUD_HANDOFF.md. The joint Gram backend exists and passes toy independent tests; no actual-source joint numerical proposal or accepted combined energy path exists yet.

No research compute process is live. The former subagents have completed; spawn new bounded agents if helpful. The source folder is not a Git repository. Existing source code and all result/history files remain intact. A local zip backup exists in `/Users/aidenlippert/Documents/Spectra-cloud-transfer-20260912/`; it was not uploaded. The unused cloud task is not the research destination.

Start by confirming the checkpoint and running the focused charged tests, then make actual progress on the coupled Gram attack. Do not stop after restating this document.
