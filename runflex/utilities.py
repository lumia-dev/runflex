#!/usr/bin/env python

import os
import subprocess
from loguru import logger
from typing import Union
from pathlib import Path
from importlib import metadata


def getfile(filename: str) -> Path:
    files = [_.locate() for _ in metadata.files('runflex') if str(_).endswith(filename)]
    assert len(files) == 1, logger.critical(f"Can't resolve the path to {filename} as several files share the name")
    return files[0]


def runcmd(cmd: str):
    logger.info(cmd)
    return subprocess.run(cmd.split())


def checkpath(path: Union[Path, str]) -> Path:
    if not os.path.isdir(path):
        os.makedirs(path)
    return Path(path)