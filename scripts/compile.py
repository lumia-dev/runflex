#!/usr/bin/env python
import sys

from runflex.rctools import rc
from runflex.logging_tools import logger
from runflex.runtools import runFlexpart
from argparse import ArgumentParser, REMAINDER
from pandas import read_csv

p = ArgumentParser()
p.add_argument('--rc', help="Main configuration file (i.e. rc-file)", required=True)
p.add_argument('--verbosity', '-v', default='DEBUG')
p.add_argument('--obs', '-o', help="Observation file (CSV-formatted)", required=True)
p.add_argument('args', nargs=REMAINDER)
args = p.parse_args(sys.argv[1:])

logger.setLevel('DEBUG')

# Load the observations
db = read_csv(args.obs, parse_dates=[0], index_col=False, infer_datetime_format=True)

# Run
fp = runFlexpart(rc(args.rc))
fp.compile()
fp.setupObs(db)
fp.distribute()