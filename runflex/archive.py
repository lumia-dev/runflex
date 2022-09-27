#!/usr/bin/env python
from typing import List
from pathlib import Path
from datetime import datetime
from loguru import logger
from pandas import DataFrame
import os
from subprocess import Popen, PIPE
from tqdm import tqdm
from h5py import File
import time


class Lock:
    def __init__(self, filename: str, timeout : float = 3600, delay : float = 5.):
        self.filename = filename
        self.timeout = timeout
        self.delay = delay
        self.lock = None

    @property
    def locked(self) -> bool:
        return self.lock is not None

    def __enter__(self):
        self.acquire()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.release()

    def __del__(self):
        if self.lock is not None :
            self.release()

    def acquire(self):
        while True:
            try :
                self.lock = File(self.filename, mode='w')
                return self
            except (OSError, BlockingIOError) as e :
                if datetime.now().timestamp() - Path(self.filename).stat().st_mtime > self.timeout:
                    self.release(warn=f'Removing expired lock file {self.filename}')
                logger.info(f"Lock file found {self.filename}, waiting {self.delay} seconds")
                time.sleep(self.delay)

    def release(self, warn : str = None):
        if warn:
            logger.warning(warn)
        self.lock.close()
        self.lock = None
        Path(self.filename).unlink()


class Rclone:
    def __init__(self, address: str, lockfile: str = None, lock_expire: int = 3600, max_attempts: int=3):
        """
        :param address: location of the rclone repo
        :param lockfile: path to the lock file.
        :param lock_expire: expiration time for the lock file (avoid blocking on locks that have been left there before)
        :param max_attempts: Max number of download attempts
        """
        self.address = address
        self.lockfile = lockfile
        self.lock_expire = lock_expire
        self.max_attempts = max_attempts
        self.filelist = self.lockfile + '.list'
        self.lock = None

    def get(self, files_list: DataFrame, dest: Path, attempts_nb: int = 0, info: str = None) -> bool:
        """
        :param files_list: list of files to be downloaded
        :param dest: folder where the files need to be put
        :param attempts_nb: number of attempts
        :param info: optional message to be added to the progress bar
        :return:
        """

        # If no address has been provided, still check if the files are there
        if not self.address:
            return self.check(list(files_list.file), dest)

        # Download files one month at a time.
        # Do not pre-check the files presence, as rclone does this better.
        for files in files_list.groupby(files_list.time.dt.month):
            month = files[1].time.iloc[0].strftime('%Y/%-m')

            # Add the current month to the progress bar text
            if info is not None :
                msg = info + '; ' + month
            else :
                msg = month

            # retrieve the files:
            self.retrieve(files[1].file.values, dest, month, info=msg)

        # Check if files are there, and repeat if not.
        while not self.check(files_list.file.values, dest):
            if attempts_nb == self.max_attempts :
                return False
            return self.get(files_list, dest, attempts_nb + 1)
        return True

    def retrieve(self, files: List[str], dest: Path, source: str = '', info: str =''):
        """
        :param files: list of files to be retrieved
        :param dest: destination directory
        :param source: source directory
        :param info: optional text added to the progress bar
        :return:
        """

        with Lock(self.lockfile) as lockf:
            with open(self.filelist, 'w') as fid:
                fid.writelines(file+'\n' for file in files)

            msg = f'Downloading meteo from {self.address}'
            if info is not None :
                msg += f' ({info})'

            # adapted from https://gist.github.com/wholtz/14b3a3479fe9c70eefb2418f091d2103
            cmd = f'rclone copy -P {os.path.join(self.address, source)} --include-from {self.filelist} {dest}'.split()
            logger.info(cmd)
            with tqdm(total=100, desc=msg, unit="%") as pbar:
                with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                    for line in proc.stdout:
                        line = line.strip()
                        if line.startswith('Transferred:') and line.endswith('%'):
                            percent = float(line.split(',')[1].split('%')[0])
                            pbar.n = percent
                            pbar.refresh()

    @staticmethod
    def check(files: List[str], dest: Path) -> bool:
        files_missing = any(not Path(dest).joinpath(f).exists() for f in files)
        return not files_missing
    #
    # def __enter__(self):
    #     self.lock = Lock(self.lockfile)
    #     self.lock.acquire()
    #
    # def __exit__(self, exc_type, exc_val, exc_tb):
    #     self.lock.release()
    #
    # def __del__(self):
    #     self.lock.release()


# class Rclone_:
#     def __init__(self, address: str, lockfile: Path = None, lock_expire: int = 3600, max_attempts: int=3, task_id = ''):
#         self.address = address
#         self.lockfile = Path(lockfile)
#         self.lock_expire = lock_expire
#         self.max_attempts = max_attempts
#         self.lock = None
#
#     def get(self, files_list: DataFrame, dest: Path, attempt_nb : int = 0, flatten: bool = True, info=None) -> bool:
#         """
#         :param files_list: list of files to download
#         :param dest: folder in which the files should be put
#         :param attempt_nb: number of attempts (for internal use)
#         :param flatten: whether the files should all end up in the same subfolder, or respect the file structure of the remote repo
#         :param info: optional message to be added to the progress bar
#         :return:
#         """
#
#         all_files_present = self.check(list(files_list.file), dest)
#
#         if not self.address or all_files_present:
#             return all_files_present
#
#         if not flatten :
#             self.retrieve(files_list.file.values, dest, info=info)
#         else:
#             for files in files_list.groupby(files_list.time.dt.month) :
#                 month = files[1].time.iloc[0].strftime('%Y/%-m')
#                 if info is not None :
#                     info += '; ' + month
#                 else :
#                     info = month
#                 self.retrieve(files[1].file.values, dest, month, info=info)
#
#         while not self.check(files_list.file.values, dest) :
#             if attempt_nb == self.max_attempts :
#                 return False
#             return self.get(files_list, dest, attempt_nb + 1)
#         return True
#
#     def retrieve(self, files: List[str], dest: Path, source: str = '', info=None):
#         with open(self.lockfile, 'a') as fid :
#             fid.writelines(file+'\n' for file in files)
#
#         msg = f'Downloading meteo from {self.address}'
#         if info is not None :
#             msg += f' ({info})'
#
#         # adapted from https://gist.github.com/wholtz/14b3a3479fe9c70eefb2418f091d2103
#         cmd = f'rclone copy -P {os.path.join(self.address, source)} --include-from {self.lockfile} {dest}'.split()
#         with tqdm(total=100, desc=msg, unit="%") as pbar:
#             with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
#                 for line in proc.stdout:
#                     line = line.strip()
#                     if line.startswith('Transferred:') and line.endswith('%'):
#                         percent = float(line.split(',')[1].split('%')[0])
#                         pbar.n = percent
#                         pbar.refresh()
#
#     @staticmethod
#     def check(files: List[str], dest: Path) -> bool:
#         files_missing = any(not Path(dest).joinpath(f).exists() for f in files)
#         return not files_missing
#
#     def __enter__(self) -> "Rclone":
#         """
#         Check for the existence of a lockfile. Normally, all "Tasks" handled by a same "Manager" share the same uid, and therefore the same lockfile. So this avoids that two tasks try to retrieve data in the same time.
#         If the lockfile is too old (default: 3600 seconds), then the lock is released.
#         """
#         # nsec = 1.
#         delay = 1.
#         while True:
#             if self.lockfile.exists():
#                 if datetime.now().timestamp() - self.lockfile.stat().st_mtime > self.lock_expire:
#                     self.release_lock(warn=f'Removing expired lock file {self.lockfile}')
#                 logger.info(f"Meteo lock file found {self.lockfile}, waiting {delay} seconds")
#                 sleep(delay)
#                 delay = min(delay + 5, 30)
#             else :
#                 self.lock = open(self.lockfile, 'a')
#                 return self
#
#         # while self.lockfile.exists():
#         #     sleep(nsec)
#         #     # Get the age of the file:
#         #     try :
#         #         if datetime.now().timestamp() - self.lockfile.stat().st_mtime > self.lock_expire:
#         #             self.release_lock(warn=f'Removing expired lock file {self.lockfile}')
#         #         logger.info(f"Meteo lock file found {self.lockfile}, waiting {nsec} seconds")
#         #         nsec = min(nsec + 5, 30)
#         #     except FileNotFoundError:
#         #         # It can happen that the file has been deleted between the "while" statement and the "os.path.getctime" one.
#         #         # Then just stay in the while loop and test again.
#         #         pass
#         # open(self.lockfile, 'a').close()
#         # return self
#
#     def __exit__(self, exception_type, exception_value, traceback):
#         self.release_lock()
#
#     def release_lock(self, warn: str = None):
#         self.lock.close()
#         if warn :
#             logger.warning(warn)
#         try :
#             self.lockfile.unlink(missing_ok=True)
#         except FileNotFoundError :
#             # If for some reason, the lock has been removed by another way, just pass
#             pass
#
