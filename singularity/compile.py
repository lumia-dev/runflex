#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from runflex import compile as runtools
from runflex.logging_tools import logger

p = ArgumentParser()
p.add_argument('--builddir', '-b', help="Directory where the FLEXPART source code will be copied", required=True)
p.add_argument('--src', '-s', help='Base source directory', required=True)
p.add_argument('--verbosity', '-v', default='INFO')
p.add_argument('--extra', help="Optional extra source paths, which overwrites the source files in the base directory.")
p.add_argument('--fc', help='Fortran compiler', default='gfortran')
p.add_argument('--dest', '-d', help="Installation path", default='/usr/bin')
args = p.parse_args(sys.argv[1:])

logger.setLevel(args.verbosity)

src = runtools.Source(args.src, extras=args.extra, machine='singularity', compiler=args.fc, builddir=args.builddir).compile(args.dest)