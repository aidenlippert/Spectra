"""A checked template operation that avoids re-normalizing every concrete product.

The learned program only proposes a template over sorted rank labels. Every new
template is compared with the original exact CAR normalizer and with literal
ladder action on its <=6-mode local algebra. Only that checked template enters
the cache. Runtime substitution is injective and order-preserving.
"""
from research.rsi_discovery_20260916.v2.algebra import oracle, validate_polynomial, ladder_identity


class CertifiedTemplates:
    def __init__(self, proposer):
        self.proposer = proposer
        self.templates = {}
        self.proposal_cache = {}
        self.checks = 0
        self.hits = 0

    def __call__(self, left, right):
        for word in (left, right):
            if len(word) > 3 or any(type(c) is not int or c not in (0, 1) or type(i) is not int or i < 0 for c, i in word):
                raise ValueError("Certified template domain: valid degree<=3 CAR words")
        modes = sorted({i for _, i in left + right})
        ranks = {i: rank for rank, i in enumerate(modes)}
        pattern = (tuple((c, ranks[i]) for c, i in left), tuple((c, ranks[i]) for c, i in right))
        if pattern not in self.templates:
            candidate = self.proposer({"left": pattern[0], "right": pattern[1], "cache": self.proposal_cache})
            validate_polynomial(candidate)
            if candidate != oracle(*pattern) or not ladder_identity(*pattern, candidate):
                raise ValueError("Proposed CAR template does not have an exact witness")
            # Own the immutable result; a later proposal cannot mutate the checked witness.
            self.templates[pattern] = tuple((w, int(c)) for w, c in candidate.items())
            self.checks += 1
        else:
            self.hits += 1
        return {tuple((c, modes[i]) for c, i in word): value for word, value in self.templates[pattern]}

    def receipt(self):
        return {"checked_templates": self.checks, "checked_template_reuses": self.hits,
                "maximum_local_algebra_modes": 6, "full_molecular_sector_enumerated": False,
                "trusted_operations": ["exact CAR normalizer", "literal local ladder witness", "increasing relabeling"],
                "scope": "Each admitted template, transported to injective increasing mode labels; not arbitrary learned Python correctness"}
