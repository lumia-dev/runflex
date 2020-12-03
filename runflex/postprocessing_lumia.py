import os
import glob
from h5py import File
from netCDF4 import Dataset, chartostring
from datetime import datetime, timedelta
from pandas import DataFrame
from copy import deepcopy
from numpy import nonzero, array, meshgrid, array_equal, unique, float16, int16, float32
from multiprocessing import Pool
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


class LumiaFootprintFile:
    def __init__(self, filename, obsdb):
        self.filename = filename
        self.obsdb = obsdb
        tmin = self.obsdb.time.min()
        self.origin = datetime(tmin.year, tmin.month, 1)
        self.global_attrs = {}
        self.empty = self.init()

    def init(self):
        empty = True
        with File(self.filename, 'a') as ds:
            if 'latitudes' in ds :
                self.lats = ds['latitudes'][:]
                self.lons = ds['longitudes'][:]
                self.tres = ds.attrs['tres']
                self.origin = datetime.strptime(ds.attrs['start'], '%Y-%m-%d %H:%M:%S')
                for attr in ds.attrs:
                    self.global_attrs[attr] = ds.attrs[attr]
                empty = False
        return empty

    def add(self, obsid, filename):
        fp = SingleFootprintFile(filename=filename)
        status = fp.check_valid()
        if status > 0 :
            return status

        # On the other hand, if anything fails in what comes next, we want the whole code to fail:
        if self.empty :
            self.lats = fp.coords['lats']
            self.lons = fp.coords['lons']
            self.tres = fp.dt.total_seconds()
            self.empty = False
            with File(self.filename, 'a') as ds :
                ds['latitudes'] = self.lats
                ds['longitudes'] = self.lons
                ds.attrs['tres'] = self.tres
                ds.attrs['start'] = self.origin.strftime('%Y-%m-%d %H:%M:%S')
        else :
            assert array_equal(self.lats, fp.coords['lats'])
            assert array_equal(self.lons, fp.coords['lons'])
            assert self.tres == fp.dt.total_seconds()

        # Make sure that all data is on the same time coordinates
        fp.shift_origin(new_origin=self.origin)

        # Store global attributes :
        release_attrs = {}
        for attr in fp.ncattrs:
            if attr.startswith('release_'):
                release_attrs[attr] = fp.ncattrs[attr]
            elif attr not in self.global_attrs :
                self.global_attrs[attr] = fp.ncattrs[attr]
            elif fp.ncattrs[attr] != self.global_attrs[attr]:
                release_attrs[attr] = fp.ncattrs[attr]

        # Write obs
        with File(self.filename, 'a') as ds :
            try :
                if obsid in ds :
                    msg = f'{ds[obsid].attrs["History"]}; Modified on {datetime.today().strftime("%Y-%m-%d %H:%M")}'
                    del ds[obsid]
                else :
                    msg = f'Created on {datetime.today().strftime("%Y-%m-%d %H:%M")}'
            except RuntimeError :
                print(self.filename)
                raise RuntimeError
            gr = ds.create_group(obsid)
            gr.attrs['History'] = msg
            gr['ilons'] = fp.data.ilon.values
            gr['ilats'] = fp.data.ilat.values
            gr['itims'] = fp.data.itim.values
            gr['sensi'] = fp.data.value.values
            for attr in release_attrs:
                gr.attrs[attr] = release_attrs[attr]
                #gr.attrs[attr] = fp.ncattrs[attr]
            for attr in self.global_attrs:
                ds.attrs[attr] = self.global_attrs[attr]
        return 0

    def concat(self, field='inpfile'):
        #print(f"Adding data to {self.filename}")
        status = []
        for obs in self.obsdb.itertuples():
            status.append(self.add(obs.obsid, getattr(obs, field)))
        self.obsdb.loc[:, 'status'] = status
        return self.obsdb.status


class LFPF(LumiaFootprintFile):
    def __init__(self, filename, origin):
        self.filename = filename
        self.origin = origin
        self.global_attrs = {}
        self.empty = self.init()


class Concat:
    def __init__(self, db, field_input='inpfile', field_output='outpfile'):
        self.db = db
        self.field_input = field_input
        self.field_output = field_output
        self.files = []
        for file in unique(db.get(field_output)):
            self.files.append(LumiaFootprintFile(file, db.loc[db.get(field_output) == file]))

    def _concat(self, fp):
        return fp.concat(field=self.field_input)

    def run(self, ncpus=None):
        p = Pool(processes=ncpus)
        status = []
        _ = [status.extend(r) for r in tqdm(p.imap(self._concat, self.files), total=len(self.files))]
        self.db.loc[:, 'status'] = status
        return array(status)


class Concat2:
    """
    Same as Concat, except that it doesn't rely on a user database:
    - the grouping is done following the site names, one file/site/month
    - the obsids are generated automatically:
    """
    def __init__(self, path, dest=None, ncpus=1, remove_files=False):
        if ncpus < 1 :
            ncpus = os.cpu_count()

        self.path = path
        dest = dest if dest is not None else path
        self.dest = dest
        self.remove_files = remove_files
        
        # Get the list of files to concatenate:
        flist = glob.glob(os.path.join(path, '*.*.??????????????.hdf'))
        self.flist = array([os.path.basename(f) for f in flist])

        # Guess their sitename and times:
        self.sitemonth = array([f"{a}.{b}.{c[:6]}" for (a,b,c,d) in [f.split('.') for f in self.flist]])
        self.obsdates = array([datetime.strptime(c.split('.')[2], '%Y%m%d%H%M%S') for c in self.flist])

        # Loop over the individual files:
        p = Pool(processes=ncpus)
        _ = [r for r in tqdm(p.imap(self.worker, unique(self.sitemonth)), total=len(unique(self.sitemonth)))]

    def worker(self, sm):

        # Create an output file with the month as origin
        site, height, month = sm.split('.')
        month = datetime.strptime(month, '%Y%m')
        outfile = os.path.join(self.dest, month.strftime(f"{site}.{height}.%Y-%m.hdf"))
        out = LFPF(outfile, month)

        # Select the list of files to add, and create their obsids:
        files = self.flist[self.sitemonth == sm]
        dates = self.obsdates[self.sitemonth == sm]
        obsid = [f"{site}.{height}.{d.strftime('%Y%m%d-%H%M%S')}" for d in dates]
        files = [os.path.join(self.path, f) for f in files]

        # Add the files to the output file:
        for (fid, oid) in zip(files, obsid):
            out.add(oid, fid)
        
        if self.remove_files :
            for fid in files :
                os.remove(fid)


class SingleFootprintFile:
    def __init__(self, **kwargs):
        self.ncattrs = {}
        self.data = DataFrame(columns=['itim', 'ilat', 'ilon', 'value'])
        self.coords = {}
        self.release_start = None
        self.release_id = None
        self.dt = None
        self.origin = None
        self.filename = None
        self.valid = False
        if 'data' in kwargs:
            self.loadData(**kwargs)
        elif 'filename' in kwargs:
            self.readHDF(kwargs.get('filename'))

    def check_valid(self):
        """
        Check the validity of the footprint file:
        - return 1 if the file couldn't even be loaded
        - return 2 if the file could be loaded but was empty
        - return 0 otherwise
        """
        if not self.valid :
            return 1
        if self.data.shape[0] == 0 :
            return 2
        else :
            return 0

    def loadData(self, ncattrs=None, release=None, coords=None, data=None, specie=None):
        self.ncattrs = ncattrs

        # Iterate over the dictionary elements, which ensures we work on a copy and not on a pointer
        for k, v in coords.items():
            self.coords[k] = deepcopy(v)

        # Load the data (keep only the non-zero in memory)
        data = data.reshape(-1)
        sel = nonzero(data)
        self.data = DataFrame({
            'itim':coords['grid'][0][sel], 
            'ilat':coords['grid'][1][sel],
            'ilon':coords['grid'][2][sel],
            'value':data[sel]}
        )

        # Set the main attributes:
        self.release_start = release.start
        self.release_id = release.id
        self.dt = self.coords['time'][1]-self.coords['time'][0]
        
        # concatenate the attribute dictionaries
        release_attrs = {
            'release_id':self.release_id,
            'release_lat':f'{(release.lat1+release.lat2)/2.:.2f}',
            'release_lon':f'{(release.lon1+release.lon2)/2.:.2f}',
            'release_height':f'{(release.z1+release.z2)/2.:.1f}',
            'release_kindz':f'{release.kindz:.0f}',
            'release_time':(release.start+(release.end-release.start)/2).strftime('%Y-%m-%d %H:%M:%S'),
            'release_npart':f'{release.npart:.0f}',
            'release_lat1':f'{release.lat1:.2f}',
            'release_lat2':f'{release.lat2:.2f}',
            'release_lon1':f'{release.lon1:.2f}',
            'release_lon2':f'{release.lon2:.2f}',
            'release_z1':f'{release.z1:.2f}',
            'release_z2':f'{release.z2:.2f}',
            'release_start':self.release_start.strftime('%Y-%m-%d %H:%M:%S'),
            'release_end':release.end.strftime('%Y-%m-%d %H:%M:%S'),
        }
         
        self.species_attr = {}
        for k, v in specie.items():
            self.ncattrs[f'specie_{k}'] = v
        for k, v in release_attrs.items():
            self.ncattrs[k] = v
        self.ncattrs['tres'] = self.dt.total_seconds()

        self.origin = self.coords['time'].min()-self.dt/2
        #self.ncattrs['origin'] = self.origin.strftime('%Y-%m-%d %H:%M:%S')

        self.valid = self.data.shape[0] > 0

    def shift_origin(self, new_origin=None):
        """
        Use the beginning of the month as convention for the origin, by default (even if that leads to negative indices)
        """
        # Set the new time origin, if it is not provided as argument
        if new_origin is None:
            new_origin = datetime(self.release_start.year, self.release_start.month, 1) #+self.dt/2

        # Determine the number of time intervals there is between the two origins
        nshift = (self.origin-new_origin)/self.dt
        assert nshift-int(nshift) == 0
        nshift = int(nshift)

        # Apply this to the self.data.itim array
        self.data.loc[:, 'itim'] = self.data.loc[:, 'itim']+nshift

        # Save the new origin:
        self.origin = new_origin

        # Construct the new time axis:
        self.coords['time'] = self.itime_to_time(range(self.data.itim.min(), self.data.itim.max()))

    def itime_to_time(self, itimes):
        return array([self.origin+(n+0.5)*self.dt for n in itimes])  # 1st time element should be half interval after the origin
        
    def writeHDF(self, path):
        fname = os.path.join(path, f'{self.release_id}.hdf')

        # Modified attributes:
        self.ncattrs['origin'] = self.origin.strftime('%Y-%m-%d %H:%M:%S')

        with File(fname, 'w') as ds :

            # netCDF attributes
            for k,v in self.ncattrs.items():
                ds.attrs[k] = v

            # Coordinates
            ds['latitudes'] = self.coords['lats'].astype(float16)
            ds['longitudes'] = self.coords['lons'].astype(int16)
            ds['times'] = array([x.timetuple()[:6] for x in self.coords['time']]).astype(int16)

            # Data:
            ds['itim'] = self.data.itim.values.astype(int16)
            ds['ilon'] = self.data.ilon.values.astype(int16)
            ds['ilat'] = self.data.ilat.values.astype(int16)
            ds['value'] = self.data.value.values.astype(float32)

    def to_gridTime(self):
        from numpy import zeros
        """
        This is a diagnostic method, written to check that we can correctly reconstruct the original "grid_time" files
        Here we completely neglect the netCDF attributes and focus on the "spec001_mr" array
        """
        # Reconstruct the original time coordinate: 
        etime = datetime.strptime(self.ncattrs['iedate']+self.ncattrs['ietime'], '%Y%m%d%H%M%S')
        btime = datetime.strptime(self.ncattrs['ibdate']+self.ncattrs['ibtime'], '%Y%m%d%H%M%S')
        nt = int((etime-btime)/self.dt)
        dt = -int(self.dt.total_seconds())

        times = range(dt, (nt+1)*dt, dt)

        # Reconstruct the FLEXPART time indices:
        times_mid = self.itime_to_time(self.data.itim.values)   # Center of the time steps
        times_beg = [t-self.dt/2 for t in times_mid]
        tt = [(t-etime).total_seconds() for t in times_beg]
        itimes = [times.index(int(t)) for t in tt]

        # Create the output array
        data_out = zeros((nt, len(self.coords['lats']), len(self.coords['lons'])))
        for it, ilat, ilon, value in zip(itimes, self.data.ilat.values, self.data.ilon.values, self.data.value.values):
            data_out[it, ilat, ilon] = value
        return {'data':data_out, 'times':times, 'lats':self.coords['lats'], 'lons':self.coords['lons']}

    def readHDF(self, fname):
        self.filename = fname
        try :
            with File(fname, 'r') as ds:

                # netCDF attributes:
                for k, v in ds.attrs.items():
                    self.ncattrs[k] = v

                # Coordinates:
                self.coords['lats'] = ds['latitudes'][:]
                self.coords['lons'] = ds['longitudes'][:]
                self.coords['time'] = array([datetime(*x) for x in ds['times']])

                # Data:
                self.data.loc[:, 'itim'] = ds['itim'][:]
                self.data.loc[:, 'ilon'] = ds['ilon'][:]
                self.data.loc[:, 'ilat'] = ds['ilat'][:]
                self.data.loc[:, 'value'] = ds['value'][:]

                self.dt = timedelta(seconds=self.ncattrs['tres'])
                self.release_start = datetime.strptime(self.ncattrs['release_start'], '%Y-%m-%d %H:%M:%S')
                self.release_id = self.ncattrs['release_id']
                if 'origin' in self.ncattrs :
                    self.origin = datetime.strptime(self.ncattrs['origin'], '%Y-%m-%d %H:%M:%S')
                else :
                    self.origin = datetime(self.release_start.year, self.release_start.month, 1)
            self.valid = True
        except OSError :
            self.valid = False


class MfpFile:
    def __init__(self, filename):
        self.filename = filename
        self.ncattrs = {}
        self.coords = {}
        self.dt = None
        self.releases = None
        self.specie = {}
        self.corrfac = 1.
        self.load_metadata()

    def load_metadata(self):
        with Dataset(self.filename) as ds :
            # Load attributes
            for attr in ds.ncattrs():
                self.ncattrs[attr] = getattr(ds, attr)

            # Load coordinates
            self.dt = timedelta(seconds=float(ds['time'][0]))
            self.t_start = datetime.strptime(self.ncattrs['ibdate']+self.ncattrs['ibtime'], '%Y%m%d%H%M%S')
            self.t_end = datetime.strptime(self.ncattrs['iedate']+self.ncattrs['ietime'], '%Y%m%d%H%M%S')
            self.coords['time'] = self._transform_time(ds['time'][:])
            self.coords['lats'] = ds['latitude'][:]
            self.coords['lons'] = ds['longitude'][:]
            surface_layer_depth = ds['height'][:]
            if len(surface_layer_depth) > 1 :
                logger.error("A one-layer footprint is expected by this script")
                raise NotImplementedError

            #self.coorfac = 28.97/1000./surface_layer_depth[0]   # m3.s/kg to m2.s/mol
            self.coorfac = ds['spec001_mr'].weightmolar/1000/surface_layer_depth[0]
            self._gen_coordinates()

            # Load release info:
            self.releases = DataFrame({
                'id': array([t.strip() for t in chartostring(ds['RELCOM'][:])]),
                'lat1':ds['RELLAT1'][:],
                'lat2':ds['RELLAT2'][:],
                'lon1':ds['RELLNG1'][:],
                'lon2':ds['RELLNG2'][:],
                'z1':ds['RELZZ1'][:],
                'z2':ds['RELZZ2'][:],
                'kindz':ds['RELKINDZ'][:],
                'start':[self.t_end+timedelta(seconds=float(tt)) for tt in ds['RELSTART'][:]],
                'end':[self.t_end+timedelta(seconds=float(tt)) for tt in ds['RELEND'][:]],
                'npart':ds['RELPART'][:],
            })

            # Load species info:
            for attr in ds['spec001_mr'].ncattrs() :
                self.specie[attr] = getattr(ds['spec001_mr'], attr)

            # Add the original filename to ncattrs, for tracability:
            self.ncattrs['filename'] = os.path.abspath(self.filename)

    def __iter__(self):
        """
        Iterate over the individual footprints and return "SingleFootprintFile" instances,
        which can then be written to a file
        """
        with Dataset(self.filename) as ds :
            # Shift the order of the time dimension, to match with the coordinates
            for irl, release in self.releases.iterrows():
                fp = SingleFootprintFile(
                    ncattrs=self.ncattrs,
                    release=release,
                    coords=self.coords,
                    data=ds['spec001_mr'][0, irl, :, 0, :, :][::-1,:,:]*self.corrfac,
                    specie=self.specie,
                )
                if fp.valid:
                    fp.shift_origin()
                yield fp

    def _transform_time(self, dtimes):
        # 1st, make sure that the simulation was backward in time
        assert self.dt.total_seconds() < 0
        
        # Convert the "dtimes" to timedelta
        dtimes = array([timedelta(seconds=float(dt)) for dt in dtimes])

        # Add to the reference time. We "remove" self.dt/2 at the end to get to the middle of the time step
        # since self.dt is negative, this advances the times array by half a time step
        times = self.t_end + dtimes - self.dt/2

        # Finally, we revert the time axis to have it in increasing order
        return times[::-1]

    def check(self, path):
        """
        Check that the reconstructed footprints are correct
        """
        with Dataset(self.filename) as ds :
            for irl, release in tqdm(self.releases.iterrows(), total=len(self.releases)):
                fp = SingleFootprintFile(filename=os.path.join(path, f'{release.id}.hdf'))
                test = fp.to_gridTime()
                if array_equal(ds['spec001_mr'][0, irl, :, 0, :, :], test['data']):
                    logger.info(f"Release {release.id} successfully reconstructed!")
                else :
                    logger.warn(f"Reconstruction of release {release.id} failed ...")
            return test

    def _gen_coordinates(self):
        """
        Just generate a grid of coordinates, based on the current coordinate arrays:
        """
        nt = len(self.coords['time'])
        nlat = len(self.coords['lats'])
        nlon = len(self.coords['lons'])
        self.coords['grid'] = [g.reshape(-1) for g in meshgrid(range(nt), range(nlat), range(nlon), indexing='ij')]


class split_gridtime:
    def __init__(self, pattern, dest, ncpus=1):
        self.dest = dest
        if not os.path.exists(dest):
            os.makedirs(dest)
        fnames = glob.glob(pattern)
        valid = [os.path.exists(os.path.join(os.path.dirname(f), 'flexpart.ok')) for f in fnames]
        fnames = array(fnames)[valid]

        if ncpus < 1 :
            ncpus = os.cpu_count()
        p = Pool(processes=ncpus)
        _ = [r for r in tqdm(p.imap(self.worker, fnames), total=len(fnames))]

    def worker(self, fname):
        fpf = MfpFile(fname)
        for fp in fpf :
            fp.writeHDF(self.dest)
        print(fname)

# def split_gridtime(filename, dest):
#     fpf = MfpFile(filename)
#     for fp in fpf:
#         fp.writeHDF(dest)


def PostProcessor(run):
    """
    This postprocessor splits the multi-footprints FLEXPART "grid_time_*.nc" files in single-footprint nc files.
    By default the new files are in the same folder as the original grid_time files. Set this to a different folder
    using the "path.output.pp" rc-key.
    """
    rcf = run.rcf

    # Get the paths
    path = rcf.get('path.output')
    dest = rcf.get('path.output.pp', default=rcf.get('path.output'))

    # Guess the output file name. There may be several "grid_time" files (older runs accidentally left there)
    # so make sure we take only the most recent
    outfiles = glob.glob(os.path.join(path, 'grid_time_*.nc'))
    outfiles.sort(key=os.path.getmtime)
    outfile = outfiles[-1]

    # For complete safety, compare the creation time of the header and grid_time files. header should be older
    # but not by more than 12 hours 
    header_age = os.path.getmtime(os.path.join(path, 'header'))
    output_age = os.path.getmtime(outfile)
    if 0 < (output_age-header_age)/3600 < 12:
        split_gridtime(outfile, dest)
    elif output_age < header_age:
        logging.warn(f"header older than grid_time file {outfile}. The simulation probably crashed: skipping post-processing")
    else :
        logging.warn(f"header suspiciously old compared to the model output in {path} ({(output_age-header_age)/3600} hours). Skipping the post-processing.")
