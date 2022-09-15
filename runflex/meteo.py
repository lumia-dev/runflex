#!/usr/bin/env python

from pandas import Timedelta, Timestamp, date_range, DataFrame
from dataclasses import dataclass
from datetime import datetime
from loguru import logger
import os
from pathlib import Path
import glob


@dataclass(kw_only=True)
class Meteo:
    path : Path
    archive : str
    tres : Timedelta
    prefix : str = 'EA'
    lockfile : Path = 'runflex.rclone.meteo.lock'

    def __post_init__(self):
        # Put the lockfile in the same directory as the meteo files, unless a specific path is provided
        if not os.path.dirname(self.lockfile):
            self.lockfile = self.path.joinpath(self.lockfile)

    def __setattr__(self, key, value):
        if key in ['tres', 'path']:
            try :
                value = self.__annotations__[key](value)
            except (TypeError, ValueError) as e:
                logger.critical(f"Can't convert value {key} (value {value}) to {type} ")
                raise e
        super().__setattr__(key, value)

    def check_unmigrate(self, start: Timestamp, end: Timestamp) -> bool:

        # Generate the list of files
        files = self.gen_filelist(start, end)

        # Attempt to clone files from archive:
        with self.archive as archive:
            return archive.get(files, self.path)

    def gen_filelist(self, start: Timestamp, end: Timestamp) -> DataFrame:
        """
        Generate a list of meteo files that FLEXPART will need.
        """
        times = date_range(start - self.tres, end + self.tres, freq=self.tres)
        return DataFrame.from_dict({'time': times, 'file': [f'{self.prefix}{tt:%y%m%d%H}' for tt in times]})

    def write_AVAILABLE(self, filepath: str) -> None:
        fmt = f'{self.prefix}????????'
        flist = glob.glob(os.path.join(self.path, fmt))
        times = [datetime.strptime(os.path.basename(f), f'{self.prefix}%y%m%d%H') for f in flist]
        with open(filepath, 'w') as fid :
            fid.writelines(['\n']*3)
            for tt in sorted(times) :
                fid.write(tt.strftime(f'%Y%m%d %H%M%S      {self.prefix}%y%m%d%H         ON DISC\n'))