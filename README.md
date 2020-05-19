# RUNFLEX

Runflex is a python library to run FLEXPART on a HPC cluster

## Getting started
git clone the source code wherever you want to have it. 

### Pre-requisites
* Python packages:
    - numpy
    - pandas
    - tqdm (optional, normally)

* For FLEXPART (see FLEXPART documentation for more details)
    - a fortran compiler
    - netcdff library
    - GRIB/ECCODES library

* Run environment:
    - interactive jobs on any bash environment
    - slurm job scheduler implemented

### Content
* the _runflex_ folder contains the runflex python library itself
* the _scripts_ folder contains **example** run script (_calcFootprints.py_) and configuration files:
    - the _options_ and _species_ folders contain the standard FLEXPART configuration files (COMMAND, SPECIES_*** and OUTGRID files)
    - the _flexpart.rc_ file is an example configuration file for runflex.
* the _flexpart10.4_ folder contains the FLEXPART10.4 source code
* the _bin_ folder contains one single "_runflex_" script, which can be used to prepare the runflex environment on tetralith (set environment variables, load the required modules, make sure that runflex is within the PYTHONPATH)

### Installation and use
- the _runflex_ folder must be within the `$PYTHONPATH` environment variable
- the actual paths, besides that of the _runflex_ library, are all set in the rc-file containing the _runflex_ settings (_config.rc_ in the example). Make sure that the `path.src` points to the right source code, that the `path.inputs` points to where the FLEXPART input data are located, etc.
- the example `config.rc` rc-file reads-in some environment variables (`$RUNFLEX_PATH`, `$INTERACTIVE`, `$RUNFLEX_NCPUS`, `$RUNFLEX_SCRATCH`) ==> it's possible to replace them by hardcoded values, but they are adapted to tetralith (and to HPCs using slurm in general), in particular they define the number of parallel runs dynamically, depending on the job settings.
- the _bin/runflex_ script automatizes setting the paths and variables on tetralith. Just call (interactively or within a job) with `runflex mycommand`, where `mycommand` can be anything, but would typically be `python`, `ipython` or a python script importing runflex (like _scripts/calcFootprints.py_). This assumes that the _bin/runflex_ script is within your bash search path (either add the _bin_ folder to your `$PATH` environment variable, or copy/link _bin/runflex_ to some place within `$PATH`)

### Example script
See the inline documentation in the example scripts