"""Use the existing canonical-word exact projector during preparation only."""
import argparse
from pathlib import Path
from research.collective_completion_20260914 import spin_rows
from research.sector_quotient_20260914.fast_twirl import twirl
from research.transfer_solver_20260915.prepare import build


def run(case):
    spin_rows.twirl = twirl
    build(case.resolve())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path)
    run(parser.parse_args().case)
