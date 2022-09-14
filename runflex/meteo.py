#!/usr/bin/env python

from pandas import Timedelta, Timestamp, date_range
from dataclasses import dataclass
from typing import List
from datetime import datetime
from loguru import logger
from time import sleep
import subprocess
import os
from pathlib import Path
import glob


@dataclass
class Archive:
    address : str
    lockfile : str = None
    lock_expire : int = 3600
    max_attempts : int = 3

    def get(self, files: List[str], dest: Path, attempt_nb : int = 0) -> bool:
        if not self.address:
            return self.check(files, dest)

        with open(self.lockfile, 'a') as fid :
            fid.writelines(file+'\n' for file in files)

        subprocess.run(['rclone', 'copy', '-P', '--files-from', self.address, self.lockfile, dest])

        while not self.check(files, dest) :
            if attempt_nb == self.max_attempts :
                return False
            return self.get(files, dest, attempt_nb + 1)
        return True

    @staticmethod
    def check(files: List[str], dest: Path) -> bool:
        files_missing = any(not os.path.exists(os.path.join(dest, f)) for f in files)
        return not files_missing

    def __enter__(self) -> "Archive":
        """
        Check for the existence of a lockfile. Normally, all "Tasks" handled by a same "Manager" share the same uid, and therefore the same lockfile. So this avoids that two tasks try to retrieve data in the same time.
        If the lockfile is too old (default: 3600 seconds), then the lock is released.
        """
        nsec = 1.
        while os.path.exists(self.lockfile):
            sleep(nsec)
            # Get the age of the file:
            try :
                if datetime.now().timestamp() - os.path.getctime(self.lockfile) > self.lock_expire:
                    self.release_lock(warn=f'Removing expired lock file {self.lockfile}')
                logger.info(f"Meteo lock file found {self.lockfile}, waiting {nsec} seconds")
            except FileNotFoundError:
                # It can happen that the file has been deleted between the "while" statement and the "os.path.getctime" one.
                # Then just stay in the while loop and test again.
                pass
            nsec = min(nsec + 5, 30)
        open(self.lockfile, 'a').close()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.release_lock()

    def release_lock(self, warn: str = None):
        if warn :
            logger.warning(warn)
        try :
            Path(self.lockfile).unlink(missing_ok=True)
            os.remove(self.lockfile)
        except FileNotFoundError :
            # If for some reason, the lock has been removed by another way, just pass
            pass


@dataclass(kw_only=True)
class Meteo:
    path : Path
    archive : str
    tres : Timedelta
    prefix : str = 'EA'
    lockfile : str = 'runflex.rclone.meteo.lock'

    def __post_init__(self):
        # Put the lockfile in the same directory as the meteo files, unless a specific path is provided
        if not os.path.dirname(self.lockfile):
            self.lockfile = os.path.join(self.path, self.lockfile)

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
        with Archive(self.archive, lockfile=self.lockfile) as archive :
            return archive.get(files, self.path)

    def gen_filelist(self, start: Timestamp, end: Timestamp) -> List[str]:
        """
        Generate a list of meteo files that FLEXPART will need.
        """
        times = date_range(start - self.tres, end + self.tres, freq=self.tres)
        return [f'{self.prefix}{tt:%y%m%d%H}' for tt in times]

    def write_AVAILABLE(self, filepath: str) -> None:
        fmt = f'{self.prefix}????????'
        flist = glob.glob(os.path.join(self.path, fmt))
        times = [datetime.strptime(os.path.basename(f), f'{self.prefix}%y%m%d%H') for f in flist]
        with open(filepath, 'w') as fid :
            fid.writelines(['\n']*3)
            for tt in sorted(times) :
                fid.write(tt.strftime(f'%Y%m%d %H%M%S      {self.prefix}%y%m%d%H         ON DISC\n'))