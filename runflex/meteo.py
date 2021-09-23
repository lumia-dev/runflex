#!/bin/env python

import os
import glob
import logging
from runflex.logging_tools import logging
from time import sleep
import subprocess
from datetime import datetime, timedelta
import shutil
from numpy import array

logger = logging.getLogger(__name__)


class LocalArchive:
    def __init__(self, key):
        self.path = key
    
    def get(self, filename, count=0):
        # Construct the full archived file path:
        dest, filename = os.path.split(filename)
        date = datetime.strptime(filename[2:], '%Y%m%d%H')
        fpath = os.path.join(self.path, str(date.year), str(date.month))
        fname = os.path.join(fpath, filename)

        # Try to get the file:
        if os.path.exists(fname):
            shutil.copy2(fname, dest)
        else :
            return False
        
        # Check if the file is now there:
        if os.path.exists(os.path.join(dest, filename)):
            return True
        elif count < 5 : # If it's still not there, retry, up to 5 times.
            sleep(count+1)
            return self.get(filename, count=count+1)
        else : # If after 5 tries, it's still not there, return False
            return False


class RcloneArchive:
    def __init__(self, key, lockfile_path='/${TMPDIR}/${USER}/', lockfile_expire=3600., lockfile='runflex.rclone.meteo.lock'):
        _, self.remote, self.path = key.split(':')
        self.lockfile = os.path.join(lockfile_path, lockfile)
        self.remote_structure = {}
        self.remote_files = []
        self.lock_expire = lockfile_expire

    def get(self, filename):
        self._acquire_lock()

        # Get the path on the archive:
        dest, filename = os.path.split(filename)
        date = datetime.strptime(filename[2:], '%y%m%d%H')
        path = os.path.join(self.path, str(date.year), str(date.month))

        # Get the list of files in the rclone repo, in the same folder (only if we are in a new folder)
        success = False
        if path not in self.remote_structure :
            logger.info(f"Retrieving list of files in rclone folder {self.remote}:{path})")
            try :
                self.remote_structure[path] = subprocess.check_output(['rclone', 'lsf', f'{self.remote}:{path}'], universal_newlines=True).split('\n')
            except subprocess.CalledProcessError :
                # Previous command will return an exception if the path doesn't exist. 
                logger.warning(f"Path {self.remote}:{path} not found on rclone archive")
                self.remote_structure[path] = []
        
        # Retrieve the file
        if filename in self.remote_structure[path] :
            cmd = ['rclone', 'copy', f'{self.remote}:{path}/{filename}', dest]
            logger.warning(' '.join(cmd))
            _ = subprocess.run(cmd)
            success = os.path.exists(dest)

        self._release_lock()
        return success

    def _acquire_lock(self):
        while os.path.exists(self.lockfile):
            # Get the age of the file:
            if datetime.now().timestamp()-os.path.getctime(self.lockfile) > self.lock_expire:
                self._release_lock(warn=f'Removing expired lock file {self.lockfile}')
            logger.warn(f"Meteo lock file found {self.lockfile}, waiting 30 seconds")
            sleep(30)
        open(self.lockfile, 'a').close()

    def _release_lock(self, warn=False):
        if warn:
            logger.warn(warn)
        try :
            os.remove(self.lockfile)
        except FileNotFoundError :
            # If for some reason, the lock has been removed by another way, just pass
            pass


class Archive:
    def __init__(self, key, *attrs, **kwattrs):
        self.key = key
        if os.path.isdir(key):
            self.type = 'local'
        elif key.startswith('rclone:'):
            self.archive = RcloneArchive(key, *attrs, **kwattrs)
        else : 
            logger.error(f"Un-recognized meteo archive: {key}")
            raise RuntimeError
        self.get = self.archive.get


class Meteo:
    def __init__(self, path, archive=None, prefix='EN', tres=None, minspace=None, minage=None, **kwargs):
        """
        path (str): folder where meteo files will be read-in by FLEXPART
        prefix (str): prefix of the meteo files (e.g. ENXXXXXXXX)
        tres (timedelta): time interval between two meteo files
        archive (str): path to a meteo archive
        """
        self.path = path
        self.prefix = prefix
        self.tres = tres
        if archive is not None :
            archive = Archive(archive, lockfile_path=path, **kwargs)
        self.archive = archive
        if not os.path.isdir(self.path):
            os.makedirs(self.path)

    def genFileList(self, ti, tf, safe=True):
        """
        Generate a list of meteo files that the FLEXPART run will look for. Arguments are:
        - ti (datetime): start date
        - tf (datetime): final date
        - safe (bool): add one day before and after the requested period, to be sure we don't miss anything (default for now)
        returns a list of meteo file filenames
        """

        # Add safety margins 
        if safe:
            ti = datetime(ti.year, ti.month, ti.day)-timedelta(1)
            tf = tf+self.tres+timedelta(1)
        fname = []
        while ti <= tf :
            fname.append(os.path.join(self.path, ti.strftime(f'{self.prefix}%y%m%d%H')))
            ti += self.tres
        return fname

    def checkUnmigrate(self, ti, tf):

        # 1) Generate the list of files
        files = self.genFileList(ti, tf)
        fails = []

        # 2) Check if the file is present
        for file in files :

            # 3) Download it from archive if it's missing
            if not os.path.exists(file):
                success = False
                if self.archive is not None :
                    success = self.archive.get(file)
                if not success :
                    fails.append(file)

        if fails :
            msg = "Not all meteo files could be retrieved:\n"
            for fail in fails:
                msg += f"   {fail}\n"
            logger.error(msg)
            raise RuntimeError

    def genAvailableFile(self, fname):
        flist = glob.glob('%s/%s*'%(self.path, self.prefix))

        # Get the times of the files on disk
        times = []
        for ff in sorted(flist):
            times.append(datetime.strptime(os.path.basename(ff), self.prefix+'%y%m%d%H'))

        # Build the AVAILABLE file, line by line:
        lines = ['\n']*3
        lines.extend([l.strftime(f'%Y%m%d %H%M%S      {self.prefix}%y%m%d%H         ON DISC\n') for l in sorted(times)])

        # Write:
        with open(fname, 'w') as fid:
            fid.writelines(lines)

    def cleanup(self, minspace, minage):
        """
        minspace (int): minimum space (in Gb) that should be on the disk. Leave to None for disabling the cleanup
        minage (int): minimum age (in hours) of the meteo files that can be deleted
        """

        minspace = minspace*1024**3

        if shutil.disk_usage(self.path).free < minspace : 
            # List all the meteo files in the folder
            files = glob.glob(os.path.join(self.path, f'{self.prefix}????????'))

            # Find their age, in hours
            ref = datetime.today().timestamp()
            ages = (array([os.path.getatime(f) for f in files])-ref)/3600.

            # Select only the files older than minage
            files = files[ages > minage]
            files = files.sort(key=os.path.getatime, reverse=False)

            # Delete files, starting from the older ones, until no files are left or until there is enough space:
            while shutil.disk_usage(self.path).free < minspace:
                os.remove(files.pop[0])
