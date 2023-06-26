# RUNFLEX

Runflex is a python library to run FLEXPART on a HPC cluster

### Content
* the _runflex_ folder contains the runflex python library itself
* the _bin_ folder contains the main executable (_runflex_)
* the _src_ folder contains the FLEXPART fortran code.
* the _inputs_ folder contains input files for FLEXPART runs (default COMMAND file, SPECIES_* files, default runflex settings (_defaults.yaml_) and some data files).
* the _docs_ folder contains documentation (also accessible through (https://lumia-dev.github.io/runflex)

### Installation

For clarity below, the _runflex_ folder here designates the folder containing the _runflex_ python library, while the _root_ folder, designates the root of the git repository (i.e. the parent of _runflex_ folder).

From the _root_ folder, just do `pip install [-e] .`. This should also install any missing python libraries (or just call `make install`).

### Use

The main script is _bin/runflex_. In the most typical case (computation of footprints for LUMIA), it is called with:
```
runflex --footprints --compile --host my_host --rc flexpart.yaml
```
`--rc` points to a configuration file in the yaml format. It sets in particular the various paths (`paths` section), the location of the observations file for which footprints should be computed, and the various FLEXPART settings (grid, species, output type, release length, etc.).

`--host` points to the name section of the yaml file to be used as a `host` section, which stores machine-specific settings and paths: this facilitates using the same yaml file on different computers.

`--compile` and `-footprints` options specifiy whether we want to compile flexpart and/or calculate footprints with it.

For more advance documentation, please have a look at https://lumia-dev.github.io/runflex/ (or in the docs folder).