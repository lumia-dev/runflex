# Usage

## Compile flexpart

There are two main ways to compile the code:
1. via `runflex --compile --rc config.yaml` (can be combined with the `--footprints` option).
2. via `make build`.

The first option is typically made for compiling FLEXPART in a directory outside the *runflex* installation path (e.g. runflex is installed in */runflex*, and you work in a project */home/mpyproj*. Use `runflex --compile` to compile flexpart in that directory). The second option is meant for deploying runflex (and FLEXPART) on a new system. Under the hood, both options use the `runflex.compile.Flexpart` class.

When compiling with `runflex` (first option above), the relevant config file settings are:

* `paths.build`: where the code should be compiled
* `paths.makefile`: which makefile to use
* `paths.src`: where is the main source code
* `paths.extras`: where is the *extra* source code: source files in that folder overwrite those in the `path.src` folder (that enables having different source codes for different projects, e.g. a *dev* code and a *stable* one).

When compiling with `make` (second option), the four paths listed above are passed via the `--build`, `--makefile`, `--src` and `--extra` arguments of the `runflex/compile.py` script (called inside the *makefile*). The default *make* behaviour is to compile two branches of the code: one with `--extra` pointing to *src/extra* (compiled inside *build/flexpart/stable*), and one with `--extra` pointing to *src/dev* (compiled inside *build/flexpart/dev*).

## Compute footprints

`runflex --footprints --rc config_file.yaml [options] --host hostname`

See the Settings section below for instruction on how to construct the *config_file.yaml* configuration file

## Singularity/Apptainer wrapper



# Settings

## YAML file and resolvers

Settings should be stored in a configuration file in *yaml* format. The file is parsed using omegaconf, which enables resolvers, i.e. the syntax `key : ${resolver:value}` will internally return `key = function(value)`, provided that the resolver "resolver" has been registered (see in runflex/config.py).

The following resolvers are pre-defined:
- `runflex`: `${runflex:value}` will return a path pointing to the file or folder "x", relative to where the "runflex" library is installed. For instance, `${runflex:inputs/COMMAND}` will point to the default "COMMAND" file.
- `rclone`: `${rclone:remote:path} will return an instance of `runflex.archive.Rclone` (`Rclone(remote, path)`)
- `tmp`: `${tmp:path}` will create a unique temporary directory inside the `path` directory.

The resolvers can be combined: `${tmp:${runflex:build}}}` will create a temporary directory under the `build` directory of the local *runflex* installation.

## --host argument

The *yaml* file should contain a `host` section pointing to the machine-specific settings (paths, etc.). However, instead of including it directly in the yaml file, it is possible to specify which section should be used as `host`, by running `runflex` with the `--host` argument:

`runflex --host laptop --rc config.yaml --footprints`

Provided that the `yaml` file has a section called `laptop`, it will be automatically renamed (in memory) as `host`. For instance, the yaml file may look like this:

```
laptop:
    paths :
        build : ${runflex:build}
        src : ${runflex:src/flexpart10.4}
        makefile : makefile.laptop
        meteo : /data/FLEXPART/meteo

hpc:
    paths :
        build : /scratch/flexpart/build
        src : ${runflex:src/flexpart10.4}
        makefile : makefile.hpc
        meteo : /proj/FLEXPART/meteo/EA5/

meteo :
    path : ${host.paths.meteo}
    ...

...
```

Depending on the value of the `--host` argument, the `host.paths.build`, `host.paths.src`, `host.paths.makefile` and `host.paths.meteo` keys will point to those defined in the `laptop` section or in the `hpc` section.

## yaml file structure

### host section

Paths specific to the local machine where flexpart is run (this particular one whould be called with the `--host laptop` argument (see above).

```
laptop :
  paths : 
    build : build
    run : run
    src : ${runflex:src/flexpart10.4}
    extras : ${runflex:src/dev}
    makefile : makefile.ubuntu.gfortran
    output : output
    logfile : ${host.paths.run}/flexpart.out
    command : ${file:COMMAND}
    meteo : /scratch/FLEXPART/meteo
```

### run section
the `serial` option can be over-written by the `--serial` argument of `runflex` (this sets ncpus to 1). 

```
run :
  logfile : ${host.paths.logfile}
  serial : no
  ncpus : 8
  releases_per_task : 50
```

### paths section
```
paths : ${host.paths}
```

The `paths` section should contain the following keys:
- run : where FLEXPART is run (can be a temporary directory)
- output : where FLEXPART final output is transferred (i.e. after postprocessing)
- command : location of the COMMAND file
- build : where FLEXPART is / should be built
- meteo : path to the meteo files

### observations
```
observations :
  coordinates :
    htm :
      start : 1 jan 2018
      end : 1 jan 2019
      freq : 3H
      lat : 56.0976
      lon : 13.4189
      alt : 115.0
      height : 150
      code : htm
      range : from 12:00 to 18:00

  release :
    kindz :
      threshold : 1000
```

The observations section can follow several conventions. Mainly, the release coordinates can be directly specified in the yaml file, as in the example above, or passed by through a file. See the `runflex.observations` module for more info.

### outgrid

How the FLEXPART "OUTGRID" file should be constructed:
```
outgrid :
  x : [-15, 35, 0.25]    # [lon_min, lon_max, lon_step]
  y : [33, 73, 0.25]     # [lat_min, lat_max, lat_step]
  levels : [100]        # Height of the surface layer
```

### releases:

Release characteristics.

The "species" sub-category should either point to a "SPECIES" file, or provide keys to construct it (using default values in the "runflex.files.SPECIES" class).

```
releases :
  length : 14
  species : 
    species : AIRTRACER
    weightmolar : 29

  mass : ${.npart}
  npart : 10_000
```

### meteo
```
meteo :
  archive : ${rclone:swestore:FLEXPART/meteo/ea.eurocom025x025}
  prefix : EA
  interv : 1h
  logfile : meteo.log
```

### postprocessing
```
postprocess :
  lumia : yes
```

### flexpart
Additional rc-keys to be passed to FLEXPART (provided that the version of FLEXPART you use supports rc-files. This is not the case with the vanilla FLEXPART10.4 branch).

```
flexpart :
  npartmax : 500_000
```