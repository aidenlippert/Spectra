# Primary-source positioning

This is a focused check, not an exhaustive priority search. No world-first
claim is supported by this pass.

* The positive ground eigenmatrix ingredient is established Hubbard
  mathematics, not a new result here. See
  [Lieb, Two Theorems on the Hubbard Model (1989)](https://doi.org/10.1103/PhysRevLett.62.1201)
  and [Boretsky, Cohn and Freericks (2017)](https://arxiv.org/abs/1712.02694).
  The latter describes the attractive-model positivity result and the
  half-filled repulsive-model spin theorem with their assumptions.
* Chebyshev filtering is an established electronic-structure eigenvalue
  technique. See the author-hosted paper
  [Zhou, Saad, Tiago and Chelikowsky (2006)](https://s2.smu.edu/yzhou/publications/scf_cheb.pdf).
  Its SCF acceleration does not imply a certified interacting many-body lower
  bound for the present model.
* Inexact Chebyshev recurrence analysis also predates this work. The
  [author's publication list](https://cims.nyu.edu/~overton/papers/itermeth.html)
  identifies Golub and Overton's 1988 convergence analysis. The finite-energy
  implication used here is derived in full in DERIVATION.md rather than
  inferred from that paper's title.
* Constructive renormalization lower bounds already exist; see
  [Kull, Schuch, Dive and Navascues (2024)](https://arxiv.org/abs/2212.03014).
  Combining local blocks or introducing another formal elimination identity
  does not establish a new method by itself.
* The broad universal problem has worst-case complexity barriers.
  [O'Gorman, Irani, Whitfield and Fefferman (2021)](https://arxiv.org/abs/2103.08215)
  prove QMA-completeness for electronic structure in a fixed basis and discuss
  Hubbard hardness on generic graphs. This does not prove that this ladder or
  realistic chemistry is hard in practice, and is not an unconditional proof
  that no general algorithm can exist.

The concrete contribution to this repository is an exact, local-tensor
implementation of a seed-overlap spectral witness, including reconstruction
of every recurrence residual. It accepts a loose interval on the same
interacting eight-site model and rejects the desired narrow interval. The
identity-operator extension beyond Hubbard is valid mathematics but has not
been implemented or shown computationally economical. Novelty and useful
scaling remain open.
