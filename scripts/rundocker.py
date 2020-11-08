#!/usr/bin/env python

from runflex.rctools import rc
from runflex.runtools import runFlexpart
import os, sys
from pandas import read_csv 
from datetime import datetime
from runflex.logging_tools import logger
from argparse import ArgumentParser

logger.setLevel("DEBUG")

p = ArgumentParser()
p.add_argument('--rc', help="Main configuration file (i.e. rc-file)", required=True)
p.add_argument('--obs', '-o', help="Observation file (CSV-formatted)", required=True)
args = p.parse_args()

# Load observations
db = read_csv(args.obs, parse_dates=[1], index_col=False)

# Load settings
rcf = rc(args.rc)

# Run 
fp = runFlexpart(rcf)
fp.setupObs(db)
fp.compile()
fp.distribute() 
