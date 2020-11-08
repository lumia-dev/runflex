#!/usr/bin/env python

import sys

try :
    from tqdm import tqdm

except :
    class Tqdm:
        def __init__(self):
            pass

        def __call__(self, iterable, *args, **kwargs):
            return iterable

        def write(self, message):
            sys.stdout.write(message)
            sys.stdout.write('\n')

    tqdm = Tqdm()
