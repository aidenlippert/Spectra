# A conditional convergence guarantee for residual targeting

Selected residual enrichment has a convergence guarantee under a certified excited-state gap. Full perturbative closure is not required for this theorem. The guarantee concerns energy convergence; it does not bound the number of retained directions, kernel coordinates, arithmetic bits, or construction time.

Let H be a finite Hermitian operator with spectrum in [b,B], W=B-b>0, ground energy E, and a certified threshold beta<=lambda_2(H), counting multiplicity. Let v be a normalized vector in U, with Rayleigh value mu<=u<beta. Decompose its residual into

\[
s=(I-P_U)Hv,\qquad t=(P_UH-\mu)v.
\]

Suppose the next space contains U and s. Write E_+ for its lowest Ritz value. Then

\[
\boxed{E_+-E\le
\left(1-\frac{\beta-u}{2W}\right)(\mu-E)
+\frac{\|t\|^2}{2W}.}
\]

If v is an exact Ritz vector, t=0. Repeating exact lowest-Ritz targeting therefore converges geometrically with factor q=1-(beta-u)/(2W)<1, as long as enrichment can continue. If a computed next vector has Rayleigh value at most E_++zeta, add zeta to the right-hand side. For errors bounded by a common a=||t||^2/(2W)+zeta, iteration gives e_k<=q^k e_0+a(1-q^k)/(1-q). This states the accuracy requirement on approximate solves explicitly.

## Proof

Put g^2=||s||^2. Orthogonality gives <v,s>=0 and <v,Hs>=g^2. The allowed trial vector v-s/W has Rayleigh decrease

\[
\mu-R(v-s/W)
=\frac{2g^2/W-(\langle s,Hs\rangle-\mu g^2)/W^2}
 {1+g^2/W^2}
\ge\frac{g^2}{2W}.
\]

Here <s,Hs>-mu g^2<=Wg^2 and g<=||(H-mu)v||<=W. Hence E_+<=mu-g^2/(2W), strictly below mu when s is nonzero.

The spectral inequality (H-E)(H-beta)>=0 holds because every excited eigenvalue is at least beta. Taking its expectation gives

\[
\|(H-\mu)v\|^2\ge(\mu-E)(\beta-\mu).
\]

Since s and t are orthogonal, g^2>= (mu-E)(beta-mu)-||t||^2. Substitution gives the claimed bound. The assumption mu<beta also rules out ground degeneracy for this threshold, since some eigenvalue lies at or below mu and at most one lies below beta.

## Checking the gap rather than assuming it

For the existing matched reference, the global complement lower c is already proved in [the size-transfer note](marginal_implicit_size_transfer.md). It alone is insufficient to bound the first excited eigenvalue: the retained reference block also contains excitations.

The new rational check computes the Jacobi LDL pivots of PH0P-cI and requires exactly one negative pivot, with no zero pivots. Inertia then proves exactly one reference-block eigenvalue below c. Combined with the complement bound, lambda_2(H0)>=c. Weyl's inequality gives beta=c-eta when ||H-H0||<=eta.

For the ten- and twelve-mode fixtures at strength 1/100, the exact inertia and witness-gap checks pass. The spectral ceiling B=binom(m,2)+m/5+eta bounds same-side repulsion plus all reference hopping terms and the perturbation. A reference ground lower minus eta supplies b. The conservative contraction factors are approximately 0.99115 and 0.99168 at ten modes, and 0.99603 and 0.99640 at twelve modes (cycle and mixed respectively). Exact parameters and rounded-outward witness uppers are recorded in `results/marginal_implicit_certificate/targeting_convergence_bound.json`.

## Relation to the implemented target chains

`targeted_workspace` takes integer coordinates for v, computes s by exact projection, and adds its entire H0 closure. It therefore contains the direction required by the proof. Since v is only an approximate Ritz vector, its internal residual t is generally nonzero. Exact projection and physical-norm rounding preserve a valid vector; they do not by themselves certify a sufficiently small internal residual or a next-vector minimization error zeta.

The exact-Ritz theorem is now connected to the implemented approximate pipeline by an exact acceptance gate. For each actual integer target a in metric G and projected Hamiltonian A, it recomputes

\[
\|t\|^2=\frac{(Aa)^T G^{-1}(Aa)}{a^TGa}-\mu^2.
\]

For the next integer vector, an exact positive-definite check of A_next-ell G_next certifies its subspace minimum exceeds ell. Hence zeta<=mu_next-ell. Short rational bounds round both errors outward before applying

\[
\frac{\|t\|^2}{2W}+\zeta\le(1-q)\epsilon_{\rm floor}.
\]

A fixed u is chosen by rounding the initial target's Rayleigh value upward; every later target and final vector must stay at or below that same u. The gap, Weyl bound, error terms, retained lower bounds, exact residual extensions, and complete target chain are reconstructed during standard-library replay. Serialized diagnostic claims are not trusted.

New targeted runs and repeated enrichment enforce these gates before exporting a certificate. The optional `targeting_guarantee` recipe records the fixed upper, error floor, and retained lower bounds. `--certify-targeting` can certify an existing unrefined chain without changing its energy endpoints. Legacy energy certificates remain valid without claiming this convergence property. A Krylov upper refinement exports its separately checked energy certificate without carrying over the retained-vector convergence recipe.

The default declared floor is 1e-7. It bounds the asymptotic error allowance of the accepted recurrence; it is not the current interval width or a claim that a capped chain has already converged. The exact recurrence gives e_k<=q^k e_0+(1-q^k)epsilon_floor. The finite energy intervals remain much sharper than this conservative global rate.

The 64-dimensional cap can stop enrichment before the conservative theorem reaches its desired accuracy, and the theorem makes no claim that a small space will suffice for general chemistry.


## Accepted gated runs

The fresh ten-mode cycle chain 6→14→20→26 has width 1.4444004615211528e-9; the ten-mode mixed chain 6→36→42→48 has width 1.3247709882217757e-9. Their construction includes the convergence gates. The two final twelve-mode chains also pass the added gates with their original energy endpoints unchanged.

The cycle M10 internal squared residuals are bounded by 1e-24 at both steps; M12 cycle bounds are 1e-24 and 3e-24. The independently certified next-vector errors are below 2.214e-10 (M10 cycle) and 9.937e-11 (M12 cycle). These fit their exact step budgets. Other recorded quantities, including mixed-case errors, are in each `targeting_convergence` receipt.

Tests reject malformed recipes, excessive internal residual, excessive next-vector error, invalid retained lower bounds, target energies above the fixed upper, zero error floors, and requests below the supported decimal proposal budget. They also check integer scale invariance and exact regeneration of an accepted saved-chain receipt.

All six new convergence exports have independent standard-library replays. The full marginal suite passes 153 tests in 49.337 seconds, including mutation of the middle target in a two-step chain.
