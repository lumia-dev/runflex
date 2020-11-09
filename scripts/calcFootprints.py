#!/usr/bin/env python

from runflex.rctools import rc
from runflex.runtools import runFlexpart
from pandas import DataFrame, date_range
from datetime import datetime
from runflex.logging_tools import logger

logger.setLevel("DEBUG")

# Create obs list
db = DataFrame(columns=['time', 'lat', 'lon', 'alt', 'height', 'code'])
db.loc[:, 'time'] = date_range(datetime(2018,1,3), datetime(2018,1,4), freq='H')
db.loc[:, 'lat'] = 56.1
db.loc[:, 'lon'] = 13.42
db.loc[:, 'alt'] = 12.
db.loc[:, 'height'] = 100.
db.loc[:, 'code'] = 'XXX'

# Parse config file
rcf = rc('config.rc')

# Initialize "runFlexpart" object, from the "runflex.runtools" module.
fp = runFlexpart(rcf)

# Setup the observations
fp.setupObs(db)

# Compile
fp.compile()

fp.distribute() 
#fp.configure()
