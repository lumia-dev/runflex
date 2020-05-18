#!/bin/env python

from numpy import array, argsort
from datetime import datetime, timedelta
import os, shutil
import time
import logging

logger = logging.getLogger(__name__)

class meteo:
    def __init__(self, path, prefix='EN'):
        self.path = path
        self.prefix = prefix
        self.checkAvailable()

    def checkAvailable(self):
        import glob, os
        flist = glob.glob('%s/%s*'%(self.path, self.prefix))
        self.files = []
        self.times = []
        for ff in sorted(flist):
            self.files.append(ff)
            self.times.append(datetime.strptime(os.path.basename(ff), self.prefix+'%y%m%d%H'))
        self.files = array(self.files)[argsort(self.times)]
        self.times = array(self.times)[argsort(self.times)]

    def OKForRun(self, ti, tf, tres):
        # check if the nearest meteo file before and after the simulation are within one time-step
        tprev = self.times[self.times<=ti].max()
        tnext = self.times[self.times>=tf].min()
        if (ti-tprev) <= tres and (tnext-tf) <= tres:
            # check if all files in between are present (no tolerance to missing files)
            tt = tprev
            present = []
            times = []
            while tt <= tnext:
                present.append(tt in self.times)
                times.append(tt)
                tt += tres
            if not False in present :
                return True
            else :
                absent = True-array(present)
                missing_times = array(times)[absent]
                for tt in missing_times:
                    logger.info('Meteo file "%s%s" missing in folder %s', self.prefix, tt.strftime('%y%m%d%H'), self.path)

    def genAvailableFile(self, path):
        self.checkAvailable()
        lines = ['\n']*3
        lines.extend([x.strftime('%Y%m%d %H%M%S      '+self.prefix+'%y%m%d%H         ON DISC\n') for x in sorted(self.times)])
        fid = open(path, 'w')
        fid.writelines(lines)
        fid.close()

    def genFileList(self, ti, tf, tres):
        t0 = datetime(ti.year, ti.month, ti.day)-timedelta(1)
        tmax = tf+tres+timedelta(1)
        fname = []
        while t0 <= tmax :
            fname.append(t0.strftime('%Y/%-m/'+self.prefix+'%y%m%d%H'))
            t0 += tres
        return fname

    def checkUnmigrate(self, ti, tf, tres, archive):
        flist = self.genFileList(ti, tf, tres)
        missing = []
        for file in flist :
            if not os.path.isfile(os.path.join(self.path, file)) :
                success = self.unMigrate(file, archive, self.path)
                if not success : missing.append(file)
        if missing :
            # Try a second time first :
            for file in missing :
                if not os.path.isfile(os.path.join(self.path, file)) :
                    success = self.unMigrate(file, archive, self.path)
                if not success :
                    raise RuntimeError('Some meteo files could not be found. Aborting!')
        # Not sure if that makes something, but since I have some performance issues, just in case, I 
        # delete whatever could accumulate and grow large ...
        del missing
        del flist

    def unMigrate(self, file, path0, path1, attempt=0):
        import hashlib
        file0 = os.path.join(path0, file)
        file1 = os.path.join(path1, os.path.basename(file))
        if os.path.isfile(file0):
            if not os.path.isfile(file1):
                #print '[MSG] Migrating file %s to %s'%(file0, path1)
                shutil.copy2(file0, file1)
                # check md5sum
            f0 = open(file0, 'rb')
            f1 = open(file1, 'rb')
            md0 = hashlib.md5(f0.read()).hexdigest()
            md1 = hashlib.md5(f1.read()).hexdigest()
            f1.close()
            attempt += 1
            if md0 != md1 and attempt < 6 :
                time.sleep(attempt)
                if attempt > 0 :
                    logger.warn('Copy of file %s failed. Retrying ...', file0)
                self.unMigrate(file, path0, path1, attempt)
            else :
                return True
        logger.error('No valid copy of meteo file %s was found', file)
        logger.error('        Directories checked:')
        logger.error('        %s', path1)
        logger.error('        %s', path0)
        return False
