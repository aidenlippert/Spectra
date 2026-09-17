# Spectra: architecture for compounding mathematical and physical discovery

Design the most ambitious technically credible architecture for Spectra to discover new mathematics, algorithms, and physical mechanisms, turn those discoveries into executable capabilities, and use them to improve subsequent discovery. The intended scientific impact is field-wide: overcoming important barriers in quantum many-body computation and transferring useful principles to other mathematical and physical problems.

Your deliverable is a deeply technical architecture and research specification. Give the central scientific mechanism, mathematical contracts, executable interfaces, evaluation design, and a buildable first implementation. Optimize for consequential new capability and cumulative research power. Architectural sophistication and scientific novelty are separate claims; justify each.

**1. Establish the real starting point.**

Inspect the current GitHub `main` checkout. Begin with `research/rsi_discovery_20260916/`, especially `v2/RESULTS.md`, `v2/DERIVATIONS.md`, `v2/records.py`, `v2/library.py`, `v2/certified.py`, `v2/workflow.py`, and the two campaign controllers. Read `migration/source_sync_20260916/REPORT.md` and its manifest to identify unavailable historical artifacts and publication limits. Follow the current reports in constructive response, constructive compression, local response, positive cones, and direct control to their implementations and checks. Locate additional evidence by its dependencies rather than reading the entire archive indiscriminately.

Produce an evidence map separating implemented behavior, independently checked results, conditional theorems, numerical evidence, failed searches, exact obstructions, missing evidence, and conjectures. Preserve accepted Spectra results. At publication, the two completed scoped RSI studies had failed their success criteria; faster inherited routines had not established improved subsequent discovery. Check for newer evidence before carrying that assessment forward. The first campaign's success condition also needs scrutiny: a positive interpretation must require successful acquisition and an evidenced contribution from the acquired capability.

**2. Choose a central discovery mechanism.**

Develop at least three substantially different architectural hypotheses. Explain what each changes about the attainable mathematics or algorithms, how it could compound, what resources it consumes, and an experiment that could falsify it. Choose one primary architecture and defend it against its strongest alternative. Support novelty comparisons with current primary literature and concrete existing systems.

Give the system freedom to invent operator languages, algebraic identities, invariants, quotient constructions, decompositions, effective theories, algorithms, proof tactics, experimental designs, and research policies. Design a principled way to extend its representation language when the existing language cannot express a useful idea. Every extension needs explicit semantics and an acceptance contract.

Consider mechanisms such as counterexample-guided theorem and program synthesis, certified changes of representation, compositional elimination, proof-guided abstraction refinement, learned decomposition, and cross-domain transport of invariants. These are candidates to compare; the architecture should follow your strongest supported mechanism.

**3. Specify the scientific objects and their semantics.**

Define a research problem as a typed object containing the mathematical or physical model, domain and assumptions, target quantity or claim, tolerances, resource envelope, and acceptance procedure. Define capabilities by their inputs, outputs, preconditions, postconditions, executable implementation, proof or evidence, cost model, dependencies, and known failure cases.

Specify a minimal intermediate representation that can express symbolic mathematics, executable algorithms, certificates, counterexamples, physical-model hypotheses, and experimental observations while preserving their different meanings. State how composition discharges assumptions and how evidence remains bound to the exact problem and environment. Supply schemas and worked instances, including a theorem used by a synthesized program and a numerical result that remains a conjecture.

Describe the evolving research state, for example

\[
R_t=(K_t,L_t,\Pi_t,G_t,V_t,B_t),
\]

where the components represent knowledge, executable capabilities, research policies, dependency structure, validation semantics, and resources. Define the actual transition operators. Explain how a discovery changes the reachable future search space or its cost. Make invalidation, repair, branching, suspension, and recovery explicit.

**4. Work through the many-body mathematical core.**

Use Spectra's central question as the first demanding domain: can the difficult correlations be represented by a manageable collection of collective structures while the remainder is eliminated or certified together?

Propose a concrete mathematical construction, with a fully stated theorem candidate and an executable verification route. Work through its assumptions, proof obligations, interfaces, error propagation, and failure modes. Address coupled systems, charge redistribution, the physical spin domain, long-range terms, induced interactions, and collective cross terms. Permit adaptive growth where the correlations require it.

If you propose recursive reduction, specify a contract such as

\[
\mathcal T_k:(H_k,\mathcal D_k,\varepsilon_k)
\mapsto(H_{k+1},W_k,\mathcal L_k,\Delta_k,C_k),
\]

and define every object, the energy or observable inequalities certified by \(C_k\), and the composition rule. Derive the constants in an error recurrence such as \(\delta_{k+1}\le a_k\delta_k+\eta_k\), or replace it with a better proved accounting rule. Explain how the required constants are computed and certified. Identify when the induced representation remains compact and what measurable structure controls its growth.

Separate a known elimination identity from a new method for constructing and verifying its useful ingredients. Charge dictionary construction, hidden sector enumeration, reference-state discovery, precision, memory, and rejected attempts. Evaluate candidate representations by accepted bounds or observables and their complete costs. Give a route for distinguishing optimizer failure from an insufficient representation, including exact family obstructions where feasible.

State the problem classes on which efficiency is claimed. Make every structural assumption testable or explicitly external; identify any verification of an assumption that is itself computationally hard. Account for bit complexity as well as arithmetic counts. Locate the proposal relative to known complexity barriers without assuming that all physically relevant instances share worst-case difficulty.

**5. Build genuine recursive improvement into the design.**

Specify how an acquired lemma, representation, algorithm, or policy causes a later search to find something it previously could not find as reliably or cheaply. Include improvement of search languages, decomposition strategies, experiment selection, and proof construction alongside ordinary implementation optimization.

Define a measurable research objective over held-out problems, with all acquisition and deployment costs charged. One possible form is

\[
J_B(R,\mathcal P)=\mathbb E_{p\sim\mathcal P}
\left[\max_{c:\,V_p(c)=\mathrm{accept},\;\mathrm{Cost}(R,p,c)\le B}u_p(c)\right],
\]

with the empty set, utility, budget, distribution, and stochastic uncertainty explicitly defined. Improve this objective if it misses enabling discoveries or long-term value. Keep operationally useful, enabling, and speculative capabilities distinct so that delayed-value ideas can survive without being counted as completed successes.

Design causal experiments across multiple generations. Include a strong frozen system, acquired-library ablations, policy ablations, matched resource envelopes, new problems, repeat acquisitions, and actual dependency or consumption evidence. Separate inherited execution speed from improved method invention using appropriate controls; explain where co-adaptation makes a common-backend comparison incomplete. Specify protections against benchmark leakage, repeated-selection bias, and changes to the scientific target. Define the statistical unit and the success criterion before confirmatory evaluation.

The decisive chain is: independently acquired capability → materially better subsequent discovery → frozen transfer → useful complete result → repetition across later generations. A useful implementation speedup remains valuable even when this stronger chain fails.

**6. Define the trusted boundary and durable implementation.**

Separate proposal generation, exploratory execution, proof construction, scientific acceptance, and experimental observation. Specify the smallest practical trusted core, the independent replay path, and the conditions under which a new checker or proof rule can be admitted. Mathematical acceptance and empirical validation require different contracts. Keep campaign acceptance rules fixed while evaluating a candidate system.

Provide interfaces and pseudocode for the complete research loop, capability installation, theorem transport, counterexample feedback, policy updates, and acceptance. Include failure semantics and an explicit example of a rejected update. Make resource bounds, process isolation, reproducibility, and dependency invalidation operational.

Design durable artifact storage with content hashes, full dependency closure, verified backups, and clean-machine recovery. Address the actual offloaded-file failures in the current worktree. A stored hash establishes identity only when its referenced contents are available and checked. Include a recovery experiment that reconstructs a retained capability from its published dependencies.

Map the architecture to existing modules: identify what to retain, deepen, replace, or retire, with a named production caller for every proposed component. Favor deep modules and narrow interfaces. Introduce new infrastructure where the scientific workflow requires it, with explicit ownership and measured costs.

**7. Make the path to field-wide breakthroughs concrete.**

Select one primary many-body bottleneck and two substantively different mathematical or physical challenges on which the central mechanism could transfer. For each, state the current best relevant methods, the unsolved barrier, the candidate new mechanism, the assumptions, the decisive experiment, and the outcome that would matter outside Spectra. Use primary sources for comparisons and novelty claims.

Explain how discoveries transfer through proved correspondences, conserved structure, reusable algorithms, or empirical hypotheses with new validation. Give at least one explicit cross-domain translation, including what fails to transfer. Select problems for scientific consequence and information gained about the mechanism; hydrogen-chain size alone is insufficient as a research objective.

Extend the architecture toward forward prediction, inverse design, and experimental planning. Specify what additional mathematical and physical capabilities each requires: identifiability, competing solutions, model discrepancy, uncertainty propagation, kinetics, manufacturability, and experimental feedback as applicable. Treat inverse design as a constrained inference and decision problem whose solution may be nonunique or nonexistent. Keep solver guarantees distinct from guarantees about real materials. Physical execution and new external spending must follow the user's actual authorization.

**8. Use independent technical challenges and return a decisive specification.**

Use parallel subagents for bounded investigations of the mathematical core, algorithm and representation search, implementation and trusted boundaries, causal evaluation, and physical significance. Give each a concrete question and evidence responsibility. Assign a skeptical reviewer to attack assumptions, hidden exponential work, novelty claims, and the proposed RSI success criterion. Resolve disagreements using proofs, counterexamples, measurements, or clearly retained uncertainty.

Return one integrated architecture, including:

1. The central thesis and the specific mechanism expected to produce compounding research capability.
2. An architecture diagram, module responsibilities, executable interfaces, and a complete discovery-to-reuse example.
3. Formal contracts, a developed theorem candidate, proof obligations, and complexity and error accounting.
4. A comparison with existing approaches that identifies precisely what is proposed as new.
5. A preregisterable multigeneration evaluation protocol and explicit negative-result interpretations.
6. Three breakthrough campaigns, each with a consequential target and a decisive falsification experiment.
7. A repository migration map and the first bounded implementation package, including files, tests, resource envelope, and acceptance criteria.
8. The strongest unresolved objection and the experiment or mathematical result needed to resolve it.

Make the architecture as ambitious as the scientific mechanism supports. Develop its difficult mathematical and algorithmic content in full. Finish with a specification from which another team could build and test the first complete research cycle, and show exactly how a successful cycle would create a more capable next one.
