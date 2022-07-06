#!/usr/bin/env python

from netCDF4 import Dataset
from numpy import unique, array
import os
from tqdm import tqdm
from loguru import logger
import tarfile
from pandas import read_csv


def check_success(df, path):
    """
    Check if the obs in the dataframe are already present in the footprint files
    """
    df = df.copy()

    # 1) Guess the names of the footprint files (to avoid going over all the files):
    filenames = df.code + '.' + df.height.astype(int).astype(str) + 'm.' + df.time.dt.strftime('%Y-%m.hdf')
    filenames = array([os.path.join(path, filename) for filename in filenames])

    # 2) Check in each file
    df.loc[:, 'present'] = False
    for filename in tqdm(unique(filenames), desc=f"Checking presence of footprints in {path}"):
        dfs = df.loc[filenames == filename]
        try :
            with Dataset(filename) as ds :
                df.loc[filenames == filename, 'present'] = [o in ds.groups for o in dfs.obsid]
        except FileNotFoundError :
            logger.warning(f'File {filename} not found')

    return df.loc[:, 'present']


def read_obsdb(fname):
    with tarfile.open(fname, 'r:gz') as tar:
        df = read_csv(tar.extractfile('observations.csv'), infer_datetime_format='%Y%m%d%H%M%S', index_col=0, parse_dates=['time'])
        if not 'code' in df.columns :
            df.loc[:, 'code'] = df.loc[:, 'site']
    return df.loc[:, ['time', 'lat', 'lon', 'alt', 'height', 'code', 'kindz', 'obsid']].drop_duplicates(subset='obsid')
