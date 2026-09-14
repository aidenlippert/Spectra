# Evaluation gates declared before the integrated run

Date: 2026-09-10. This protocol was written after component debugging, before the first integrated multi-seed result. It is not a blinded external preregistration.

Fixed evaluation seeds: 7, 11, 23. Three stages per seed, two Hamiltonian coefficient contexts per stage. Public vocabulary: every nonidentity Pauli string of weight at most two on three qubits (36 terms); the identity offset is fixed to zero. Each later stage introduces an absent interaction. Acquisition and evaluation use disjoint observable weights. All methods receive the same derivative-oracle observations and bounded noise assumption.

The pipeline gate requires all 18 contexts to recover the supported terms, maximum held-out derivative error <=0.006, exact interval consistency and an independently accepted rational energy certificate, numeric ground-energy bracketing, and a reduced-dynamics cross-check <= the analytic residual bound plus 1e-10 numerical tolerance. Non-identifiable probes, a planted weight-three interaction, and an undersized closure budget must cause abstention or rejection. A failing gate is not silently removed by changing seeds or thresholds.

The intelligence gate is stricter and independent: repeated discoveries must improve total cost or experimental efficiency over matched conventional and retrieval/cache baselines on new coupled families. Count acquisition, candidate search, fit, verification, and reuse. Lower correlation-search counts alone do not pass this gate. Neither a generic fixed dictionary, scripted experiment bank, nor algebraic closure computed after support recovery demonstrates a newly acquired research operation. This implementation is expected to provide a baseline and may produce a null intelligence result.

The exact certificates are conditional mathematical statements. This experiment does not provide a shot-based apparatus model, finite-difference derivative error theorem, unknown vocabulary discovery, a polynomial many-body solver, material fabrication or a physical experiment. These missing capabilities must remain visible in the result.
