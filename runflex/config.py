#!/usr/bin/env python

from omegaconf import OmegaConf, DictConfig
from runflex import prefix
from runflex.utilities import getfile
from runflex.archive import Rclone
from runflex.manager import uid
from pandas import Timedelta


OmegaConf.register_new_resolver("prefix", lambda x: prefix.joinpath(x))
OmegaConf.register_new_resolver("file", lambda x: getfile(x))
OmegaConf.register_new_resolver("rclone", lambda x: Rclone(x, lockfile=f'rclone.meteo.{uid}.lock'))
OmegaConf.register_new_resolver('dt', lambda x: Timedelta(x))
