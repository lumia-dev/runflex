#!/usr/bin/env python

import os
import subprocess
from loguru import logger
from typing import Union, List
from pathlib import Path
from importlib import metadata


#def getfile(filename: str) -> Path:
#    files = [_.locate() for _ in metadata.files('runflex') if str(_).endswith(filename)]
#    assert len(files) == 1, logger.critical(f"Can't resolve the path to {filename} as several files share the name")
#    return Path(files[0])

def getfile(filename: str) -> Path:
    files = [_.locate() for _ in metadata.files('runflex') if str(_).endswith(filename)]
    if len(files) == 1:
        return Path(files[0])
    # Fast fallback: check known locations (no globbing)
    repo_candidate = Path(__file__).resolve().parent.parent / "inputs" / filename
    if repo_candidate.exists():
        return repo_candidate

    share_candidate = Path(sys.prefix) / "share" / "flexpart" / "inputs" / filename
    if share_candidate.exists():
        return share_candidate

    # Keep the original error style
    logger.critical(
        f"Can't resolve the path to {filename}. "
        f"metadata matches={len(files)}, repo={repo_candidate}, share={share_candidate}"
    )
    raise FileNotFoundError(filename)


def runcmd(cmd: Union[str, List[str]]):
    logger.info(cmd)
    if isinstance(cmd, str):
        cmd = cmd.split()
    return subprocess.run(cmd)


def checkpath(path: Union[Path, str]) -> Path:
    if not os.path.isdir(path):
        os.makedirs(path)
    return Path(path)
