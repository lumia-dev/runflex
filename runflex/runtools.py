#!/usr/bin/env python

import os
import shutil
import time
import logging
from pandas import DataFrame
from numpy import inf
from datetime import datetime, timedelta
import runflex.logging_tools
import subprocess
from runflex import tqdm
import tempfile
import runflex.rctools as rctools
import runflex.meteo as meteo

logger = logging.getLogger(__name__)

def checkpath(path):
    if not os.path.isdir(path):
        os.makedirs(path)


class Namelists:
    def __init__(self, filename=None, names=None):
        self.lists = []
        if filename is not None :
            for name in names :
                self.lists = Namelist(name=name, file=filename)
    
    def addList(self, nmlist):
        self.lists.append(nmlist)

    def write(self, filename):
        for nmlist in self.lists :
            nmlist.write(filename, mode='a')


class Namelist:
    def __init__(self, name=None, file=None):
        self.keys = {}
        self.name = name 
        if file is not None :
            self.parse(file, name)
    
    def parse(self, filename, name=None):
        # Read all the lines in the file

        if name is not None : 
            self.name = name
        assert self.name is not None, logger.error("A name for the namelist must be provided")

        with open(filename, 'r') as fid :
            end = False
            while not end :
                line = fid.readline()
                if line.startswith(f'&{self.name}'): # We found the start of the namelist
                    line = fid.readline()
                    while not line.startswith(' /'): 
                        key, val = line.strip().strip(',').split('=')
                        self.keys[key.strip()] = val.strip()
                        line = fid.readline()
                    end = True
        
    def add(self, key, value, fmt=None):
        if fmt == str :
            self.keys[key] = f'"{value}"'
        elif fmt is None :
            self.keys[key] = value
        else :
            self.keys[key] = f'{value:{fmt}}'

    def write(self, filename, mode='w'):
        with open(filename, mode) as fid :
            fid.write(f'&{self.name}\n')
            for key in self.keys :
                fid.write(f' {key} = {self.keys[key]} ,\n')
            fid.write(' /\n\n')


class Command:
    def __init__(self, rcf, start, end):
        self.rcf = rcf
        self.start = start
        self.end = end
    
    def genFiles(self):
        rundir = self.rcf.get('path.run')
        self.gen_pathnames(rundir)
        self.gen_COMMAND(rundir)
        self.gen_OUTGRID(rundir)
        self.copyFiles(rundir)
        self.gen_SPECIES(rundir)

    def copyFiles(self, rundir):
        shutil.copy(self.rcf.get('landuse.file'), rundir)
        shutil.copy(self.rcf.get('landuse.z0.file'), rundir)
        shutil.copy(self.rcf.get('surfdepo.file'), rundir)

    def gen_pathnames(self, rundir):
        with open(os.path.join(rundir, 'pathnames'), 'w') as fid:
            fid.write('%s/\n' % rundir) 
            fid.write('%s/\n' % self.rcf.get('path.output'))
            fid.write('%s/\n' % self.rcf.get('path.meteo'))
            fid.write(os.path.join(rundir, 'AVAILABLE'))
            fid.write('\n')
            fid.write('=====\n')
            fid.write('=====\n')
            fid.write('=====\n')

    def gen_COMMAND(self, rundir):
        command = Namelist(file=self.rcf.get('file.command'), name='COMMAND')
        command.add('IBDATE', self.start.strftime("%Y%m%d"))
        command.add('IBTIME', self.start.strftime("%H%M%S"))
        command.add('IEDATE', self.end.strftime("%Y%m%d"))
        command.add('IETIME', self.end.strftime("%H%M%S"))
        command.write(os.path.join(rundir, 'COMMAND'))

    def gen_OUTGRID(self, rundir):
        gridfile = self.rcf.get('file.grid')
        shutil.copy(gridfile, rundir)

    def gen_SPECIES(self, rundir):
        checkpath(os.path.join(rundir, 'SPECIES'))
        nspec = self.rcf.get('species.number')
        specdir = self.rcf.get('path.species')
        for ispec in range(nspec):
            spec_index = self.rcf.get(f'species.{ispec+1}.index')
            shutil.copy(os.path.join(specdir, f'SPECIES_{spec_index:03.0f}'), os.path.join(rundir, 'SPECIES'))


class Observations:
    def __init__(self, rcf):
        self.rcf = rcf
        self.observations = Observations

    def setup(self, obslist):
        """
        Store the observation dataFrame in an attribute, and create an index
        The dataFrame must have columns 'lon', 'lat', 'alt', 'time', 'code' and 'height'
        The index is created as 'code.height.time' and is used to define the footprint file names
        """
        self.observations = obslist
        self.observations.loc[:, 'name'] = ['%s.%im.%s'%(o.code.lower(),o.height,o.time.strftime('%Y%m%d%H%M%S')) for o in self.observations.itertuples()]
        self.observations.set_index('name', inplace=True)
        if not 'kindz' in self.observations.columns :
            self.observations.loc[:, 'kindz'] = self.rcf.get('releases.kindz')
        if not 'release_height' in self.observations.columns :
            # Plain sites
            self.observations.loc[self.observations.kindz == 1, 'release_height'] = self.observations.loc[:, 'height']

            # Mountain top sites
            alt_corr = self.rcf.get('releases.altitude_correction', default=1)
            sampling_height = self.observations.loc[self.observations.kindz == 2, 'height'].values
            sampling_alt = self.observations.loc[self.observations.kindz == 2, 'alt'].values
            release_alt = sampling_height + alt_corr*(sampling_alt-sampling_height)
            self.observations.loc[self.observations.kindz == 2, 'release_height'] = release_alt
        
        self.observations.loc[:, 'mass'] = self.rcf.get('releases.mass')
        self.observations.loc[:, 'npart'] = self.rcf.get('releases.npart')

    def gen_RELEASES(self, rundir):
        # Header
        nspec     = self.rcf.get('species.number')
        specs     = []
        for ispec in range(nspec):
            specs.append(self.rcf.get(f'species.{ispec+1}.index'))

        releases = Namelists()
        ctrl = Namelist('RELEASES_CTRL')
        ctrl.add('NSPEC', nspec)
        ctrl.add('SPECNUM_REL', ", ".join([str(s) for s in specs]))
        releases.addList(ctrl)

        for obs in self.observations.itertuples():
            rl = Namelist('RELEASE')
            rl.add('IDATE1', obs.time.strftime('%Y%m%d'))
            rl.add('ITIME1', obs.time.strftime("%H%M%S"))
            rl.add('IDATE2', obs.time.strftime('%Y%m%d'))
            rl.add('ITIME2', obs.time.strftime("%H%M%S"))
            rl.add('LON1', obs.lon)
            rl.add('LON2', obs.lon)
            rl.add('LAT1', obs.lat)
            rl.add('LAT2', obs.lat)
            rl.add('Z1', obs.release_height)
            rl.add('Z2', obs.release_height)
            rl.add('ZKIND', obs.kindz)
            rl.add('MASS', obs.mass)
            rl.add('PARTS', obs.npart)
            rl.add('COMMENT', obs.Index, fmt=str)
            releases.addList(rl)
        releases.write(os.path.join(rundir, 'RELEASES'))

    def write(self, path, ncpus=1, nobsmax=None, maxdt=7):
        """
        For multi-process jobs, we want to distribute the footprint computations over many FLEXPART runs.
        This method splits the database in chunks and writes it to temporary files
        """

        # If we don't want to split
        if nobsmax is None :
            fname = os.path.join(path, 'observations.hdf')
            self.observations.to_hdf(fname)
            return [fname]

        # 1st, make sure the observations are sorted by time (so that jobs will contain as much as possible data in similart time interval)
        self.observations.sort_values('time', inplace=True)

        # Compute the number of chunks, based on the number of observations and obs/run:
        nobstot = self.observations.shape[0]
        nchunks = nobstot/nobsmax + (nobstot%nobsmax > 0)

        # If there are more CPUs than observation chunks, reduce the number of obs/chunk:
        if nchunks < ncpus :
            nchunks = ncpus
            nobsmax = nobstot/nchunks + (nobstot%nchunks > 0)

        logger.debug("    Number of CPUs detected: %i", ncpus)
        logger.debug("     Number of observations: %i", nobstot)
        logger.debug("    Max number of obs/chunk: %i", nobsmax)

        dbfiles = []
        i0 = 0
        pbar = tqdm(total=self.observations.shape[0], desc='splitting observation database')
        checkpath(self.rcf.get('path.scratch.global'))
        while i0 < self.observations.shape[0]:

            # select the slice of the observation database for that chunk
            # Make sure that the time interval between the first and last obs is not too long
            db = self.observations.iloc[i0:i0+nobsmax, :]
            db = db.loc[db.time <= db.iloc[0].time + timedelta(maxdt),:]

            # Write the database slice to a file accessible by all nodes:
            fid, fname = tempfile.mkstemp(dir=path, prefix = 'observations.', suffix = '.hdf')

            # Open and close the file descriptor to avoid "too many open files" error
            os.fdopen(fid).close()

            db.to_hdf(fname, 'obsdb')
            dbfiles.append(fname)
            i0 += db.shape[0]
            pbar.update(i0)
        pbar.close()

        return dbfiles

    def __getattr__(self, item):
        return self.observations.loc[:, item]


class runFlexpart:
    def __init__(self, rcf):
        self.rcf = rcf
        self.observations = Observations(self.rcf)

    def setupObs(self, obslist):
        self.observations.setup(obslist)

    def config_times(self):
        """ Time boundaries of the simulation 
        """
        tmax = self.observations.time.max()
        tmin = self.observations.time.min()
        lenmax = self.rcf.get('releases.length')
        tmin = tmin-timedelta(days=lenmax)
        return tmin, tmax

    def config_meteo(self, start, end):
        """ Check the meteorological files """
        prefix = self.rcf.get('meteo.prefix')
        checkpath(self.rcf.get('path.meteo'))
        mm = meteo.meteo(self.rcf.get('path.meteo'), prefix)
        tres = self.rcf.get('meteo.interv')
        tres = timedelta(tres / 24.)
        if self.rcf.haskey('meteo.archive'):
            mm.checkUnmigrate(start, end, tres, self.rcf.get('meteo.archive'))
        mm.genAvailableFile('%s/AVAILABLE' % self.rcf.get('path.run'))

    def compile(self):
        t0 = datetime.now()

        # Copy files to build path
        builddir = self.rcf.get('path.build')
        checkpath(builddir)
        srcdir = os.path.normpath(self.rcf.get('path.src')+'/')
        os.system('rsync -avh %s/*.f90 %s'%(srcdir,builddir))
        shutil.copy('%s/makefile.%s.%s'%(srcdir,self.rcf.get('machine'),self.rcf.get('compiler')),'%s/Makefile'%builddir)

        # make
        prevdir = os.getcwd()
        os.chdir(builddir)
        if os.path.exists('flexpart.x'):
            os.remove('flexpart.x')
        os.system('make')
        t1 = datetime.fromtimestamp(os.path.getmtime('flexpart.x'))
        os.chdir(prevdir)
        if t0 > t1 :
            logger.error("Compilation aborted. The whole universe is against you. Stop working and go get a drink. Seriously!")

    def configure(self):
        # Create the run directory
        rundir = self.rcf.get('path.run')
        builddir = self.rcf.get('path.build')
        checkpath(rundir)
        checkpath(builddir)

        # Copy the executable to the run directory
        shutil.copy(os.path.join(builddir, 'flexpart.x'), rundir)

        # Setup the meteo files
        start, end = self.config_times()
        #self.config_meteo(start, end)

        c = Command(self.rcf, start, end)
        c.genFiles()
        self.observations.gen_RELEASES(rundir)

    def distribute(self, nobsmax=None, maxdt=7):
        # Determine the size of the flexpart runs (nobsmax can be taken taken from rc-file, provided as argument,
        # or otherwise a default value of 50 obs/run is taken
        if nobsmax is None :
            nobsmax = self.rcf.get('nreleases.max', default=50)

        dbfiles = self.observations.write(
            path=self.rcf.get('path.scratch.global'),
            nobsmax=nobsmax, 
            ncpus=self.rcf.get('ntasks_parallell'),
            maxdt=maxdt
        )
        self.runTasks(dbfiles)

    def runTasks(self, obsfiles, chunks=None):
        """
        Submit the individual FLEXPART runs to SLURM.
        To avoid saturating the queue, tasks are submitted only when there is at least one free CPU

        The "chunks" optional argument can be used to run only some of the tasks (in debug or restart context)
        """

        if os.path.exists(os.path.join(self.rcf.get('path.run'), 'flexpart.ok')):
            os.remove(os.path.join(self.rcf.get('path.run'), 'flexpart.ok'))
        if self.rcf.get('run.interactive'):
            self.submit = self.submit_interactive
        else :
            self.submit = self.submit_sbatch

        pids = []
        for ichunk, dbf in enumerate(obsfiles) :
            # Count the number of active tasks
            if chunks is None or ichunk in chunks :
                pids.append(self.submit(dbf, ichunk))
                time.sleep(5)

        # Wait for the runs to finish
        if not self.rcf.get('run.interactive'):
            [pid.wait() for pid in pids]

    def submit_sbatch(self, dbf, ichunk):
        logpath = self.rcf.get('path.logs')
        checkpath(logpath)
        rundir = os.path.join(self.rcf.get('path.run'), '%i'%ichunk)
        cmd = [
            'srun',
            '--exclusive', '-N', '1', '-n', '1',
            '--job-name=flexpart.%i' % ichunk,
            '-o', '%s/job.%i.out' % (logpath, ichunk),
            '-e', '%s/job.%i.out' % (logpath, ichunk),
            'python', os.path.abspath(__file__), '--db', dbf, '--rc', os.path.join(self.rcf.dirname, self.rcf.filename), '--tag', self.rcf.get('tag'), '--path', rundir
        ]
        logger.info(' '.join([x for x in cmd]))

        # delay the submission if there are too many tasks running:
        ntasks = subprocess.check_output(['squeue', '-s', '-j', os.environ['SLURM_JOBID']]).count(os.environ['SLURM_JOBID'])
        ncpus = self.rcf.get('ntasks.parallell')

        logger.debug("Running tasks: %i", ntasks)
        while ntasks >= ncpus :
            logger.debug("Too many tasks (%i). Waiting 1 minute ...", ntasks)
            time.sleep(60)
            ntasks = subprocess.check_output(['squeue', '-s', '-j', os.environ['SLURM_JOBID']]).count(os.environ['SLURM_JOBID'])

        return subprocess.Popen(cmd)

    def submit_interactive(self, dbf, ichunk):
        rundir = os.path.join(self.rcf.get('path.run'), '%i'%ichunk)
        cmd = 'python %s --db %s --rc %s --path %s'%(os.path.abspath(__file__), dbf, os.path.join(self.rcf.dirname, self.rcf.filename), rundir)
        logger.info(cmd)
        os.system(cmd)

    def run(self):
        os.chdir(self.rcf.get('path.run'))
        subprocess.check_call(['./flexpart.x'])

def parse_options(args):
    from optparse import OptionParser, OptionGroup
    p = OptionParser(usage='%prog [options] rc_file')
    p.add_option('-r', '--rc', dest='rcf')
    p.add_option('-d', '--db', dest='obs')
    p.add_option('-p', '--path', dest='rundir')
    p.add_option('-t', '--tag', dest='tag')
    (options, outargs) = p.parse_args(args)
    return (options, outargs)


if __name__ == '__main__' :
    """
    For SLURM runs on multiple nodes, a separate python process must be launched on each node to access the local resources (local hard drive)
    For the sake of simplicity, this module can act as it's own executable
    
    usage : python runtools.py rcfile obsfile rundir
    obsfile is a hdf5-saved pandas dataFrame. rundir overwrites the "path.run" variable in rcfile.
    """

    from pandas import read_hdf
    import sys
    # parse arguments
    
    logger.info(os.getcwd())
    opts, args = parse_options(sys.argv)

    rcfile = opts.rcf
    dbfile = opts.obs
    rundir = opts.rundir

    # adjust the run path
    rcf = rctools.rc(rcfile)
    rcf.setkey('path.run', rundir)
    rcf.setkey('tag', opts.tag)

    # Read the obs database
    db = read_hdf(dbfile)

    # Initialize the flexpart run
    fp = runFlexpart(rcf)
    fp.setupObs(db)
    fp.configure()

    # Run, collect (i.e. gather the individual footprints in monthly files) and delete temp files
    fp.run()
