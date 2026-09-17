"""New-algorithm discovery prompts are generated from operational library manifests."""
import json

from research.rsi_discovery_20260916.v2.algebra import BASE_SOURCE


def prompt(tasks, available_methods, action, previous=None, feedback=None):
    text = """Discover a reusable exact algorithm for directly constructing a symmetric SOS Gram coefficient map.
This is a new algorithm task, not tuning numerical precision or a search-grid parameter.
Input payload['words'] is a list of arbitrary words of degree <=3, each a tuple of (creation_flag, mode).
The dictionary has homogeneous particle-number charge; order may be arbitrary. payload['cache'] is empty.
Output is dict canonical_word -> nonempty dict (i,j): nonzero Python int, with 0<=i<=j<len(words).
For diagonal coordinate (i,i), emit the canonical polynomial B_i^dagger B_i.
For i<j, emit B_i^dagger B_j + B_j^dagger B_i. Exact cancellations must be removed.
Creation operators precede annihilation operators; each group has ascending modes with fermionic signs.
The empty word represents identity. This map is used by the real non-enumerating Spectra SOS constructor.
Available base operation: oracle(left,right) -> exact dict canonical_word: rational integer coefficient.
The oracle has a 100000 concrete-pair LRU cache, cleared before each cold run.
You may rediscover/implement arbitrary pure integer algebra, direct structural factorization, transpose
identities, relabeling or compilation. Do not assume a particular support graph, mode count or dictionary ordering.
Existing baseline already includes spin-twirl caching and paired anticommutator contraction elsewhere;
neither is a new discovery claim for this task. This task is the general unpaired symmetric Gram map.
Objective: reduce complete cold construction time while preserving every exact coefficient and cross term.
Input generation, activation and independent checking are also charged. Cached answers keyed to tasks are forbidden.
All four experiment variants receive the same target, base model, development data and hard resource limits.
The available library below is selected automatically by the method registry. If empty, you may still
derive those algorithms yourself using the base operations. A transfer test follows after the method is frozen.
"""
    text += "\nDevelopment task specifications:\n" + json.dumps(tasks)
    text += "\nApplicable operational methods:\n" + json.dumps(available_methods)
    text += "\nResearch action selected for this live attempt: " + action
    text += "\nCurrent program:\n" + (previous or BASE_SOURCE)
    if feedback:
        text += "\nAll development observations available to this attempt:\n" + json.dumps(feedback)
    return text


def primitive_manifest(proposal, validation):
    return {"name": "car_kernel", "call": "car_kernel({'left': left, 'right': right, 'cache': payload['cache']})",
            "inputs": "two arbitrary degree<=3 CAR words and an initially empty shared mutable cache",
            "outputs": "dict canonical_word: nonzero integer", "implementation": proposal["source"],
            "domain": "general degree-six CAR algebra; injective mode relabelings",
            "assumptions": proposal["assumptions"], "failure_cases": proposal["failure_cases"],
            "evidence": validation, "cost_regime": "cache creation charged; tested on collective support dictionaries",
            "usage": "This callable is installed automatically in the candidate environment."}
