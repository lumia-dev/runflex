#!/usr/bin/env python

from omegaconf import OmegaConf
from runflex import prefix
from runflex.utilities import getfile
from runflex.archive import Rclone
from pandas import Timedelta
import tempfile


OmegaConf.register_new_resolver("prefix", lambda x: prefix.joinpath(x))
OmegaConf.register_new_resolver("file", lambda x: getfile(x))
OmegaConf.register_new_resolver("rclone", lambda x: Rclone(x))
OmegaConf.register_new_resolver('dt', lambda x: Timedelta(x))
OmegaConf.register_new_resolver('tmp', lambda x: tempfile.TemporaryDirectory(prefix='runflex', dir=x).name)
