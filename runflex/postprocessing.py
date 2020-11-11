import os
from h5py import File
from netCDF4 import Dataset, chartostring
from datetime import datetime, timedelta
from pandas import DataFrame
from copy import deepcopy
from numpy import nonzero, array, meshgrid
import logging
import glob

logger = logging.getLogger(__name__)

class Footprint:
    def __init__(self, ncattrs, release, coords, data, specie):
        self.ncattrs = ncattrs
        self.coords = coords
        self.release = release
        data = data.reshape(-1)
        sel = nonzero(data)
        self.data = DataFrame({
            'itim':coords['grid'][0][sel], 
            'ilat':coords['grid'][1][sel],
            'ilon':coords['grid'][2][sel],
            'value':data[sel]}
        )
        
        self.release_attrs = {
            'release_id':self.release.id,
            'release_lat':f'{(self.release.lat1+self.release.lat2)/2.:.2f}',
            'release_lon':f'{(self.release.lon1+self.release.lon2)/2.:.2f}',
            'release_height':f'{(self.release.z1+self.release.z2)/2.:.1f}',
            'release_kindz':f'{self.release.kindz:.0f}',
            'release_time':(self.release.start+(self.release.end-self.release.start)/2).strftime('%Y-%m-%d %H:%M:%S'),
            'release_npart':f'{self.release.npart:.0f}',
            'release_lat1':f'{self.release.lat1:.2f}',
            'release_lat2':f'{self.release.lat2:.2f}',
            'release_lon1':f'{self.release.lon1:.2f}',
            'release_lon2':f'{self.release.lon2:.2f}',
            'release_z1':f'{self.release.z1:.2f}',
            'release_z2':f'{self.release.z2:.2f}',
            'release_start':self.release.start.strftime('%Y-%m-%d %H:%M:%S'),
            'release_end':self.release.end.strftime('%Y-%m-%d %H:%M:%S')
        }
         
        self.species_attr = {}
        for k, v in specie.items():
            self.ncattrs[f'specie_{k}'] = v
        for k, v in self.release_attrs.items():
            self.ncattrs[k] = v
        
    def shift_origin(self, new_t0=None):
        """
        Use the beginning of the month as convention for the origin, by default (even if that leads to negative indices)
        """
        dt = self.coords['time'][1]-self.coords['time'][0]
        cur_t0 = self.coords['time'].min()
        if new_t0 is None:
            new_t0 = datetime(self.release.start.year, self.release.start.month, 1)+dt/2
        nshift = (cur_t0-new_t0).total_seconds()/dt.total_seconds()
        assert nshift-int(nshift) == 0
        nshift = int(nshift)
        self.data.loc[:, 'itim'] = self.data.loc[:, 'itim']+nshift
        newtaxis = []
        tt = new_t0
        while tt <= self.coords['time'].max():
            newtaxis.append(tt)
            tt+=dt
        self.coords['time'] = array(newtaxis)
        
    def writeHDF(self, path):
        fname = os.path.join(path, f'{self.release.id}.hdf')
        with File(fname, 'w') as ds :

            # netCDF attributes
            for k,v in self.ncattrs.items():
                ds.attrs[k] = v

            # Coordinates
            ds['latitudes'] = self.coords['lats']
            ds['longitudes'] = self.coords['lons']
            ds['times'] = [x.timetuple()[:6] for x in self.coords['time']]

            # Data:
            ds['itim'] = self.data.itim.values
            ds['ilon'] = self.data.ilon.values
            ds['ilat'] = self.data.ilat.values
            ds['value'] = self.data.value.values


class FlexpartOutput:
    def __init__(self, filename):
        self.ncattrs = {}
        self.footprints = {}
        with Dataset(filename) as ds :
            for attr in ds.ncattrs():
                self.ncattrs[attr] = getattr(ds, attr)
                
            # Convert time from model units (i.e. seconds from the simulation end) to datetime
            dt = timedelta(seconds=float(ds['time'][0]))
            if dt.total_seconds() < 0 :
                ref = datetime.strptime(self.ncattrs['iedate']+self.ncattrs['ietime'], '%Y%m%d%H%M%S')
            else :
                print("The footprint is not backward ...")
                raise NotImplementedError
                
            # Generate the coordinate axis
            times = array([ref+timedelta(seconds=float(tt)) for tt in ds['time'][:]])-dt/2
            lats = ds['latitude'][:]
            lons = ds['longitude'][:]
            
            # Reverse the time axis:
            times = times[::-1]
            
            # Generate a grid of coordinates
            grid = meshgrid(range(len(times)), range(len(lats)), range(len(lons)), indexing='ij')
            
            # Save all the coordinates in a dictionary
            coords = {
                'time': times,
                'lons': lons,
                'lats': lats,
                'grid': [g.reshape(-1) for g in grid]
            }
            
            # Store release and species info in dictionaries
            releases = DataFrame({
                'id': array([t.strip().replace("'","") for t in chartostring(ds['RELCOM'][:,1:])]),
                'lat1':ds['RELLAT1'][:],
                'lat2':ds['RELLAT2'][:],
                'lon1':ds['RELLNG1'][:],
                'lon2':ds['RELLNG2'][:],
                'z1':ds['RELZZ1'][:],
                'z2':ds['RELZZ2'][:],
                'kindz':ds['RELKINDZ'][:],
                'start':[ref+timedelta(seconds=float(tt)) for tt in ds['RELSTART'][:]],
                'end':[ref+timedelta(seconds=float(tt)) for tt in ds['RELEND'][:]],
                'npart':ds['RELPART'][:],
            })

            specie = {}
            for attr in ds['spec001_mr'].ncattrs() :
                specie[attr] = getattr(ds['spec001_mr'], attr)

            # Loop over the individual footprints
            for irl, release in releases.iterrows():
                self.footprints[release[0]] = Footprint(
                    self.ncattrs,
                    release,
                    deepcopy(coords),
                    ds['spec001_mr'][0, irl, ::-1, 0, :, :],
                    specie,
                )

    def writeHDF(self, dest):
        for fp in self.footprints.values():
            fp.shift_origin()
            fp.writeHDF(dest)

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
    header_age = os.path.getmtime(os.path.join('path', 'header'))
    output_age = os.path.getmtime(outfile)
    if 0 < (output_age-header_age)/3600 < 12:
        FlexpartOutput(output_age).writeHDF(dest)
    elif output_age < header_age:
        logging.warn(f"header older than grid_time file {outfile}. The simulation probably crashed: skipping post-processing")
    else :
        logging.warn(f"header suspiciously old compared to the model output in {path} ({(output_age-header_age)/3600} hours). Skipping the post-processing.")
        