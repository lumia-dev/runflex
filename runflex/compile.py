import os
import shutil
import subprocess
from datetime import datetime
from runflex.runtools import checkpath, safecopy, Namelist, Namelists
from runflex import meteo
from datetime import timedelta
import logging
import tempfile

logger = logging.getLogger(__name__)

def symlink(source, dest):
    try :
        os.symlink(source, dest)
    except FileExistsError :
        if os.path.isdir(dest):
            symlink(source, os.path.join(dest, os.path.basename(source)))
        elif os.path.islink(dest):
            if os.readlink(dest) != source :
                raise FileExistsError
        else :
            raise FileExistsError


class Source:
    def __init__(self, src, extras=None, builddir='/flexpart', machine=os.environ['machine'], compiler='gfortran'):
        self.src = src
        self.extras = extras
        self.builddir = builddir
        self.machine = machine
        self.compiler = compiler

    def compile(self, dest):
        t0 = datetime.now()

        # Copy files to the build path:
        checkpath(self.builddir)
        os.system(f'rsync -avh {self.src}/*.f90 {self.builddir}')
        shutil.copy(f'{self.src}/makefile.{self.machine}.{self.compiler}',f'{self.builddir}/Makefile')

        # Additional set of source files (typically, those that are specific to a machine or project, and not on the git repository)
        if self.extras is not None :
            os.system(f'rsync -avh {self.extras}/*.f90 {self.builddir}')
            try :
                shutil.copy(f'{self.extras}/makefile.{self.machine}.{self.compiler}', f'{self.builddir}/Makefile')
            except FileNotFoundError :
                pass

        # Make
        prevdir = os.getcwd()
        os.chdir(self.builddir)
        if os.path.exists('flexpart.x'):
            os.remove('flexpart.x')
        os.system('make')
        t1 = datetime.fromtimestamp(os.path.getmtime('flexpart.x'))
        if dest is not None and dest != self.builddir :
            if os.path.exists(os.path.join(dest, 'flexpart.x')):
                os.remove(os.path.join(dest, 'flexpart.x'))
            shutil.move('flexpart.x', dest)
        os.chdir(prevdir)
        if t0 > t1 :
            logger.error("Compilation aborted. The whole universe is against you. Stop working and go get a drink. Seriously!")


class Task:
    def __init__(self, rcf):
        self.end = None
        self.start = None
        self.command = None
        self.releases = None
        self.grid = None
        self.obs = None
        self.src = None
        self.species = []
        self.rcf = rcf
        self.rundir = self.rcf.get('path.run')
        self.uid = next(tempfile._get_candidate_names())

    def setup(self, start, end):
        self.setupFolders()
        self.setupExecutable()
        self.setupMeteo(start, end)
        self.gen_COMMAND(start, end)
        self.gen_OUTGRID()
        self.gen_SPECIES()
        self.gen_pathnames()
        self.gen_RELEASES()

    def setupFolders(self):
        """
        Create the input and output directories, with the flexpart executable. Compile if needed ...
        """
        # Create the run and output directories
        checkpath(self.rundir)
        checkpath(self.rcf.get('path.output'))
        symlink(self.rcf.get('landuse.file'), self.rundir)
        symlink(self.rcf.get('landuse.z0.file'), self.rundir)
        symlink(self.rcf.get('surfdepo.file'), self.rundir)

    #def setupExecutable(self, compile=False):
    #    # Compile or copy the executable
    #    if not self.rcf.haskey('flexpart.executable') or compile:
    #        self.src = Source(
    #            self.rcf.get('path.src.base'),
    #            self.rcf.get('path.src.extras', default=None),
    #            self.rcf.get('path.build', default='/flexpart'),
    #            self.rcf.get('machine', default=os.environ['machine']),
    #            self.rcf.get('compiler', default='gfortran'))
    #        self.src.compile()
    #    safecopy(self.rcf.get('flexpart.executable'), self.rcf.get('path.run'))

    def gen_COMMAND(self, start, end, **kwargs):
        self.command = Namelist(file=self.rcf.get('file.command'), name='COMMAND')
        # Override values using the optional arguments:
        for k, v in kwargs.items():
            self.command.add(k, v)

        # Fix the times so that we have conforming time steps
        dt = timedelta(seconds=int(self.command.keys['LOUTAVER']))
        start, end = start, end
        if dt < timedelta(days=1):
            start = datetime(start.year, start.month, start.day)
            end = datetime(end.year, end.month, end.day)
            while start + dt < start:
                start += dt
            while end < end:
                end += dt
        else:
            logger.error("LOUTAVER longer than 24 hours is not implemented in runflex (but it should be doable)")
            raise NotImplementedError
        self.start = start
        self.end = end

        self.command.add('IBDATE', self.start.strftime("%Y%m%d"))
        self.command.add('IBTIME', self.start.strftime("%H%M%S"))
        self.command.add('IEDATE', self.end.strftime("%Y%m%d"))
        self.command.add('IETIME', self.end.strftime("%H%M%S"))
        self.command.write(os.path.join(self.rundir, 'COMMAND'))

    def gen_OUTGRID(self, **kwargs):
        self.grid = Namelist(name='OUTGRID')
        x0, x1, dx = self.rcf.get('outgrid.x')
        y0, y1, dy = self.rcf.get('outgrid.y')
        nx = (x1-x0)/dx
        ny = (y1-y0)/dy
        self.grid.add('OUTLON0', float(x0))
        self.grid.add('OUTLAT0', float(y0))
        self.grid.add('NUMXGRID', int(nx))
        self.grid.add('NUMYGRID', int(ny))
        self.grid.add('DXOUT', float(dx))
        self.grid.add('DYOUT', float(dy))
        self.grid.add('OUTHEIGHTS', self.rcf.get('outgrid.levels'))
        self.grid.write(os.path.join(self.rundir, 'OUTGRID'))

    def gen_SPECIES(self):
        checkpath(os.path.join(self.rundir, 'SPECIES'))
        nspec = self.rcf.get('species.number')
        specdir = self.rcf.get('path.species')
        for ispec in range(nspec):
            spec_index = self.rcf.get(f'species.{ispec+1}.index')
            shutil.copy(os.path.join(specdir, f'SPECIES_{spec_index:03.0f}'), os.path.join(self.rundir, 'SPECIES'))

    def gen_pathnames(self):
        with open(os.path.join(self.rundir, 'pathnames'), 'w') as fid:
            fid.write(f'{self.rundir}\n')
            fid.write(f'{self.rcf.get("path.output")}/\n')
            fid.write(f'{self.rcf.get("path.meteo")}/\n')
            fid.write(os.path.join(self.rundir, 'AVAILABLE'))

    def gen_RELEASES(self):
        self.releases = Namelists()
        ctrl = Namelist('RELEASES_CTRL')
        ctrl.add('NSPEC', len(self.species))
        ctrl.add('SPECNUM_REL', ", ".join([str(s) for s in self.species]))
        self.releases.addList(ctrl)

        for obs in self.obs.itertuples():
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
            rl.add('ZKIND', int(obs.kindz))
            rl.add('MASS', obs.mass)
            rl.add('PARTS', int(obs.npart))
            rl.add('COMMENT', obs.Index, fmt=str)
            self.releases.addList(rl)

        self.releases.write(os.path.join(self.rundir, 'RELEASES'), mode='w')

    def compile(self):
        raise NotImplementedError

    def setupMeteo(self, start, end):
        """ Check the meteorological files """
        if not self.rcf.haskey('meteo.lockfile'):
            self.rcf.setkey('meteo.lockfile', f'runfex.rclone.meteo.lock.{self.uid}')

        self.meteo = meteo.Meteo(
            self.rcf.get('path.meteo'),
            archive=self.rcf.get('meteo.archive', default=None),
            prefix = self.rcf.get('meteo.prefix'),
            tres=timedelta(self.rcf.get('meteo.interv')/24.),
            lockfile=self.rcf.get('meteo.lockfile', default='runfex.rclone.meteo.lock')
        )

        if self.rcf.get('meteo.cleanup', default=False):
            self.meteo.cleanup(
                self.rcf.get('meteo.cleanup.minspace'),
                self.rcf.get('meteo.cleanup.minage')
            )

        self.meteo.checkUnmigrate(start, end)
        self.meteo.genAvailableFile('%s/AVAILABLE' % self.rcf.get('path.run'))

    def setupObservations(self, obslist, obsid='obsid'):
        """
        Store the observation dataFrame in an attribute, and create an index
        The dataFrame must have columns 'lon', 'lat', 'alt', 'time', 'code' and 'height'
        The index is created as 'code.height.time' and is used to define the footprint file names
        """
        self.obs = obslist
        if obsid not in self.obs.columns:
            self.obs.loc[:, obsid] = [f'{o.code.lower()}.{o.height:.0f}.{o.time.strftime("%Y%m%d%H%M%S")}' for o in self.obs.itertuples()]
        self.obs.set_index(obsid, inplace=True)

        if "kindz" not in self.obs.columns:
            self.obs.loc[:, 'kindz'] = self.rcf.get("releases.kindz")

        if "release_height" not in self.obs.columns:
            # Plain sites
            self.obs.loc[self.obs.kindz == 1, 'release_height'] = self.obs.loc[self.obs.kindz == 1, 'height']

            # Mountain top sites
            alt_corr = self.rcf.get('releases.altitude_correction', default=1)
            sampling_height = self.obs.loc[self.obs.kindz == 2, 'height'].values
            sampling_alt = self.obs.loc[self.obs.kindz == 2, 'alt'].values
            release_alt = sampling_height + alt_corr * (sampling_alt - sampling_height)
            self.obs.loc[self.obs.kindz == 2, 'release_height'] = release_alt

        self.obs.loc[:, 'mass'] = self.rcf.get('releases.mass')
        self.obs.loc[:, 'npart'] = self.rcf.get('releases.npart')

        # Species:
        for ispec in range(self.rcf.get('species.number')):
            self.species.append(self.rcf.get(f'species.{ispec + 1}.index'))

    def run(self):
        os.chdir(self.rcf.get('path.run'))
        subprocess.check_call('flexpart.x')


if __name__ == '__main__':

    import sys
    from argparse import ArgumentParser

    p = ArgumentParser()
    p.add_argument('--src', '-s', default='/runflex/flexpart10.4')
    p.add_argument('--extras', '-e', default='/runflex/extras')
    p.add_argument('--builddir', '-b', default='/flexpart')
    p.add_argument('--machine', '-m', default='singularity')
    p.add_argument('--compiler', '-c', default='gfortran')
    options = p.parse_args(sys.argv[1:])

    src = Source(options.src,
                 extras=options.extras,
                 builddir=options.builddir,
                 machine=options.machine,
                 compiler=options.compiler
                 )

    src.compile()