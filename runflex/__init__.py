#!/usr/bin/env python

import sys

try :
    from tqdm.autonotebook import tqdm
except :
    class tqdm:
        def tqdm(self, iterable, *args, **kwargs):
            return iterable
        def write(self, message):
            sys.stdout.write(message)
            sys.stdout.write('\n')
