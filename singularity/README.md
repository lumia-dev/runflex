# Singularity environment for FLEXPART

This folder contains the code to build and deploy runflex in a singularity container:
- The `singularity.def` definition files contains the recepie for building a container containing a runflex setup, i.e. with FLEXPART (compiled), the `runflex` python library installed, a set of pre-defined scripts (essentially for computing footprints), and some useful/essential tools (rclone, curl, task-spooler, etc.)
- The `runflex` python executable is a wrapper to execute commands in the container (it essentially handles mounting local directories in the container).
- The `genconfig.py` script creates a `runflex.ini` file under `${HOME}/.config/`, which contains the (machine-specific) default mount paths (for `/meteo`, `/output` and `/scratch`).

The container has interfaces to perform computations with runflex (e.g., produce footprints for LUMIA), but it also includes its own source code, which means that the entire runflex source code can be distributed as a single container image file (SIF file).

## Installation
### From the source code

1. Edit (if needed) the default paths in the `Makefile`
2. Build the container with `make build` (you need to have administrator rights to do this)
3. Generate the `runflex.ini` file with `make config`
4. Copy the `runflex` script to your user `bin` directory (default `${HOME}/.local/bin`) with `make install`

### From a singularity image file

If you only need to run the container (e.g., you want to produce footprints on a HPC, but won't do actual code development on that HPC), you can simply install the container with the following command:
`/path/to/runflex.sif install \
    --scratch /path/to/scratch \
    --meteo /path/to/meteo \
    --output /path/to/output \
    --bin /path/to/bin`

If you intend to develop the code (or just look at it), you can extract it with `/path/to/runflex.sif extract runflexdir` (with `runflexdir` the path where the code will be copied). You can then edit the code, re-compile the container and re-install it following the instructions from the previous section.

## Usage