#!/usr/bin/env python

from pandas import DataFrame
from typing import List
from runflex.files import Release, ReleasesHeader
from loguru import logger


class Releases(DataFrame):
    _metadata = ['species']
    species: List[str] = None

    @property
    def releases(self) -> Release:
        for obs in self.itertuples():
            yield Release(
                IDATE1=obs.time_start, ITIME1=obs.time_start,
                IDATE2=obs.time_end, ITIME2=obs.time_end,
                LAT1=obs.lat, LON1=obs.lon, Z1=obs.release_height,
                ZKIND=obs.kindz, MASS=obs.mass, PARTS=obs.npart,
                COMMENT=obs.obsid #Kristian april 2024: add ZZ and zweight here so they are added to the release files upon compiling
            )

    def write(self, filename: str):
        ReleasesHeader(len(self.species), self.species).write(filename, mode='w', name='RELEASES_CTRL')
        for release in self.releases:
            release.write(filename, mode='a', name='RELEASE')
