#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from runflex import compile as runtools
from runflex.rctools import rc
from datetime import datetime
from runflex.logging_tools import logger

p = ArgumentParser()
p.add_argument('--verbosity', '-v', default='INFO')
p.add_argument('--start', '-s', help="Start time", required=True)
p.add_argument('--end', '-e', help="Final time", required=True)
args = p.parse_args(sys.argv[1:])

logger.setLevel(args.verbosity)

rcf = rc()

start = datetime.strptime(args.start, '%Y%m%d')
end = datetime.strptime(args.end, '%Y%m%d')

# Default keys in singularity container
rcf.setkey('path.run', '/scratch/')
rcf.setkey('landuse.file', '/runflex/data/IGBP_int1.dat')
rcf.setkey('landuse.z0.file', '/runflex/data/surfdata.t')
rcf.setkey('surfdepo.file', '/runflex/data/surfdepo.t')
rcf.setkey('flexpart.executable', '/flexpart/flexpart.x')
rcf.setkey('path.meteo', '/meteo')
rcf.setkey('path.output', '/scratch/run')
rcf.setkey('path.src.base', '/runflex/flexpart10.4')
rcf.setkey('path.src.extras', '/runflex/extras')
rcf.setkey('machine', 'singularity')
rcf.setkey('flexpart.executable', '/usr/bin/flexpart.x')

# Additional keys:
rcf.setkey('meteo.archive', 'rclone:swestore:FLEXPART/meteo/${meteogrid}')
rcf.setkey('meteogrid', 'ea.eurocom025x025')
rcf.setkey('meteo.prefix', 'EA')
rcf.setkey('meteo.interv', 1)
rcf.setkey('file.command', '/runflex/scripts/options/COMMAND')

run = runtools.Task(rcf)
run.setupFolders()
#run.setupExecutable(compile=False)
run.setupMeteo(start, end)
run.gen_COMMAND(start, end, IOUT=0)
run.gen_pathnames()
run.run()