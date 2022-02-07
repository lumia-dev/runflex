#!/usr/bin/env python

import sys
import configparser
from argparse import ArgumentParser
import os

p = ArgumentParser()
p.add_argument('--container', default='.')
p.add_argument('--exec', default='~/.local/bin')
p.add_argument('--scratch', default=None)
p.add_argument('--output', default=None)
p.add_argument('--meteo', default=None)
p.add_argument('--extras', action='append')
args = p.parse_args(sys.argv[1:])

c = configparser.ConfigParser()
c['runflex'] = {
    'DefaultContainer': os.path.abspath(args.container),
    'scratch': args.scratch,
    'meteo': args.meteo,
    'output': args.output
}

if args.extras is not None:
    c['extras'] = {}
    for extra in args.extras :
        external, internal = extra.split(':')
        c['extras'][external] = internal

with open(os.path.join(os.environ['HOME'], '.config/runflex.ini'), 'w') as configfile:
    c.write(configfile)