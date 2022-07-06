#!/usr/bin/env python

"""
This script is aimed at computing a single flexpart footprint (for test runs). It compiles its own flexpart code, so it's meant for use in a dev container.
"""

import shutil
import sys
from runflex import utilities
from runflex.rctools import rc
from runflex.runtools import runFlexpart
from pandas import DataFrame, Timestamp
from runflex.logging_tools import logger
from argparse import ArgumentParser
from runflex.utilities import check_success
from runflex import compile as runtools


p = ArgumentParser()
p.add_argument('--verbosity', '-v', default='INFO')
p.add_argument('--rc', help="Main configuration file (i.e. rc-file)", required=True)
p.add_argument('--command', '-c', help='COMMAND file (optional). Alternatively, it can be specified as a rc-file option (file.command). The file must be in a path accessible by the singularity container')
p.add_argument('--setkey', action='append', help="use to override some rc-keys")
p.add_argument('--cleanup', action='store_true', default=False, help="Ensure that the rundir is clear from previous runs (set to False by default as this will erase anything in the scratch dir, even if it doesn't belong to runflex!)")
p.add_argument('--continue', action='store_true', default=True, help="Check the output folder to see if footprints are already present. Process only those who aren't. Set to False to avoid this behaviour", dest='continue_')
p.add_argument('--lat', type=float)
p.add_argument('--lon', type=float)
p.add_argument('--alt', type=float)
p.add_argument('--time')
p.add_argument('--height', type=float)
p.add_argument('--length')
p.add_argument('--code')
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
rcf.setkey('ntasks.parallell', 1)
rcf.setkey('collect_fail', True)
rcf.setkey('path.fail', '/output/failed/')

src = runtools.Source('/runflex/flexpart10.4', extras='/runflex/extras', machine='singularity', compiler='gfortran', builddir='/flexpart').compile('/flexpart')

if args.command :
    rcf.setkey('path.command', args.command)

rcf.setkey('run.interactive', True)

if args.setkey :
    for kv in args.setkey :
        k, v = kv.split(':')
        rcf.setkey(k, v)

if args.cleanup :
    shutil.rmtree(rcf.get('path.run'), ignore_errors=True)
    shutil.rmtree(rcf.get('path.output'), ignore_errors=True)

coords = {'lon':args.lon, 'lat':args.lat, 'alt':args.alt, 'time':Timestamp(args.time), 'code':args.code, 'height':args.height}
db = DataFrame(columns=coords.keys())
db.loc[0, :] = coords

fp = runFlexpart(rcf)
fp.setupObs(db)
fp.configure()
fp.run()