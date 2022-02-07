#!/usr/bin/env python
import os
import shutil
import sys
from runflex.rctools import rc
from runflex.runtools import runFlexpart
from pandas import read_csv
from runflex.logging_tools import logger
from argparse import ArgumentParser
from runflex.postprocessing_lumia import Concat2

p = ArgumentParser()
p.add_argument('--verbosity', '-v', default='INFO')
p.add_argument('--rc', help="Main configuration file (i.e. rc-file)", required=True)
p.add_argument('--obs', help="Observation file (csv format)", required=True)
p.add_argument('--interactive', '-i', action='store_true', default=False)
p.add_argument('--command', '-c', help='COMMAND file (optional). Alternatively, it can be specified as a rc-file option (file.command). The file must be in a path accessible by the singularity container')
p.add_argument('--ncpus', '-n', help='Number of parallell processes', default=1, type=int)
p.add_argument('--only', action='append', help="run only this site (add several times the argument for several sites)")
p.add_argument('--nobs', default=None, help="Use this to limit the number of observations (i.e. for test purposes)", type=int)
p.add_argument('--setkey', action='append', help="use to override some rc-keys")
p.add_argument('--cleanup', action='store_true', default=False, help="Ensure that the rundir is clear from previous runs (set to False by default as this will erase anything in the scratch dir, even if it doesn't belong to runflex!)")
args = p.parse_args(sys.argv[1:])

logger.setLevel(args.verbosity)

rcf = rc(args.rc)

# In a singularity container, we can hardcode some values:
rcf.setkey('path.scratch_global', '/scratch')
rcf.setkey('path.run', '/scratch/run')
rcf.setkey('path.output', '/scratch/run')
rcf.setkey('path.output.pp', '/output')
rcf.setkey('path.meteo', '/meteo')
rcf.setkey('run.tsp', 'T')
rcf.setkey('path.build', '/flexpart')
rcf.setkey('landuse.file', '/runflex/data/IGBP_int1.dat')
rcf.setkey('landuse.z0.file', '/runflex/data/surfdata.t')
rcf.setkey('surfdepo.file', '/runflex/data/surfdepo.t')
rcf.setkey('path.species', '/runflex/scripts/species')
rcf.setkey('ntasks.parallell', args.ncpus)
rcf.setkey('collect_fail', True)
rcf.setkey('path.fail', '/output/failed/')

if args.command :
    rcf.setkey('path.command', args.command)

if args.interactive :
    rcf.setkey('run.interactive', True)

if args.setkey :
    for kv in args.setkey :
        k, v = kv.split(':')
        rcf.setkey(k, v)

if args.cleanup :
    shutil.rmtree(rcf.get('path.run'), ignore_errors=True)
    shutil.rmtree(rcf.get('path.output'), ignore_errors=True)


db = read_csv(args.obs, parse_dates=[0], index_col=False, infer_datetime_format=True)
if args.only :
    select = [o in args.only for o in db.code]
    db = db.loc[select]

if args.nobs is not None :
    if db.shape[0] > args.nobs :
        db = db.iloc[:args.nobs]

fp = runFlexpart(rcf)
fp.setupObs(db)
fp.distribute()

# Concatenate output files:
if rcf.get('postprocess.lumia'):
    Concat2(rcf.get('path.output.pp'), ncpus=args.ncpus, remove_files=rcf.get('pp.remove_files', default=True))