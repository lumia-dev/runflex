#!/usr/bin/env python

from typing import List
from pathlib import Path
from datetime import datetime
from loguru import logger
from time import sleep
#from runflex.utilities import runcmd
from pandas import DataFrame
import os
from subprocess import Popen, PIPE
from tqdm import tqdm


class Rclone:
    def __init__(self, address: str, lockfile: Path = None, lock_expire: int = 3600, max_attempts: int=3):
        self.address = address
        self.lockfile = Path(lockfile)
        self.lock_expire = lock_expire
        self.max_attempts = max_attempts

    def get(self, files_list: DataFrame, dest: Path, attempt_nb : int = 0, flatten: bool = True) -> bool:
        if not self.address:
            return self.check(files_list, dest)

        if not flatten :
            self.retrieve(files_list.file.values, dest)
        else:
            for files in files_list.groupby(files_list.time.dt.month) :
                self.retrieve(files[1].file.values, dest, files[1].time.iloc[0].strftime('%Y/%-m'))

        while not self.check(files_list.file.values, dest) :
            if attempt_nb == self.max_attempts :
                return False
            return self.get(files_list, dest, attempt_nb + 1)
        return True

    def retrieve(self, files: List[str], dest: Path, source: str = ''):
        with open(self.lockfile, 'a') as fid :
            fid.writelines(file+'\n' for file in files)

        # adapted from https://gist.github.com/wholtz/14b3a3479fe9c70eefb2418f091d2103
        cmd = f'rclone copy -P {os.path.join(self.address, source)} --include-from {self.lockfile} {dest}'.split()
        with tqdm(total=100, desc=f"Downloading meteo from {self.address}", unit="%") as pbar:
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    line = line.strip()
                    if line.startswith('Transferred:') and line.endswith('%'):
                        percent = float(line.split(',')[1].split('%')[0])
                        pbar.n = percent
                        pbar.refresh()

        #runcmd(f'rclone copy -P {os.path.join(self.address, source)} --include-from {self.lockfile} {dest}')

    @staticmethod
    def check(files: List[str], dest: Path) -> bool:
        files_missing = any(not Path(dest).joinpath(f).exists() for f in files)
        return not files_missing

    def __enter__(self) -> "Rclone":
        """
        Check for the existence of a lockfile. Normally, all "Tasks" handled by a same "Manager" share the same uid, and therefore the same lockfile. So this avoids that two tasks try to retrieve data in the same time.
        If the lockfile is too old (default: 3600 seconds), then the lock is released.
        """
        nsec = 1.
        while self.lockfile.exists():
            sleep(nsec)
            # Get the age of the file:
            try :
                if datetime.now().timestamp() - self.lockfile.stat().st_mtime > self.lock_expire:
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
            self.lockfile.unlink(missing_ok=True)
        except FileNotFoundError :
            # If for some reason, the lock has been removed by another way, just pass
            pass

