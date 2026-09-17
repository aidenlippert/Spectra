"""Full replay with the experimentally tested exact paired-contraction expander."""
import argparse
import json
from pathlib import Path
import sys
from experiments import marginal_symbolic
from research.acceptance_channels_20260915.paired_exact import expand_squares,CALLS
from research.acceptance_channels_20260915.campaign import OUT


def run(case,tag,name):
    if Path(tag).name!=tag or Path(name).name!=name or tag in ('.','..') or name in ('.','..'):
        raise ValueError('Local proposal and new output names required')
    marginal_symbolic.expand_squares=expand_squares
    from research.acceptance_channels_20260915.replay_candidate import run as complete
    complete(tag,case,name)
    forbidden=[key for key in ('numpy','scipy','flint','numba','cupy','quimb','pyscf','cvxpy') if key in sys.modules]
    if forbidden:raise ValueError(('Non-standard-library accepting import',forbidden))
    record={'backend':'Python_standard_library_exact_paired_quartic_contractions',
        'all_other_accepting_steps_unchanged':True,'calls':CALLS,'forbidden_imports':forbidden,
        'same_molecular_receipt_comparison_with_original_required':True}
    with (OUT/'replays'/name/'paired_expansion.json').open('x') as stream:json.dump(record,stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('tag');p.add_argument('name')
    a=p.parse_args();run(a.case.resolve(),a.tag,a.name)
