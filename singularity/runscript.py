#!/usr/bin/env python3
import sys
import os
import subprocess
import shutil
from argparse import ArgumentParser
import configparser
import logging

logger = logging.getLogger(__name__)

p = ArgumentParser(add_help=False)
authorized_commands = ['python3', 'ipython3', 'bash']
p.add_argument('action', choices=['help', 'footprints', 'extract', 'install', 'blh', 'compile']+authorized_commands)
p.add_argument('--bin', default=os.path.join(os.environ['HOME'], '.local/bin'))
p.add_argument('--meteo', default=None, help='path where the meteo files will be downloaded')
p.add_argument('--scratch', default=None, help='path where temporary files can be written')
p.add_argument('--output', default=None, help="path where the model output will be written")
p.add_argument('--dest', default=os.environ['SINGULARITY_CONTAINER'], help='location of the container (use only during the installation)')
args, remainder = p.parse_known_args(sys.argv[1:])

if args.action == 'help':
    p.print_help()

if args.action == 'footprints':
    subprocess.run(['python3', '/runflex/singularity/calcfootprints.py']+remainder)

elif args.action == 'blh':
    subprocess.run(['python3', '/runflex/scripts/blh.py']+remainder)

elif args.action == 'compile':
    subprocess.run(['python3', '/runflex/singularity/compile.py', '-b', '/flexpart/build', '--src', '/runflex/flexpart10.4', '--extra', '/runflex/extras/', '--fc', 'gfortran', '--dest', '/flexpart'])

elif args.action == 'extract':
    shutil.copytree('/runflex', remainder[0], dirs_exist_ok=True)

elif args.action in authorized_commands :
    subprocess.run([args.action]+remainder)
