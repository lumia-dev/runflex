#!/usr/bin/env python

from omegaconf import OmegaConf, DictConfig
from runflex import prefix
from runflex.utilities import getfile


OmegaConf.register_new_resolver("prefix", lambda x: prefix.joinpath(x))
OmegaConf.register_new_resolver("pkg", lambda x: getfile(x))