# plume-tracker

Workflow for detecting plumes in CitcomS global mantle-flow models.

# Prerequisites

## Linux
`scons`: A python build utility can be installed on linux systems using the system package manager, e.g. `apt` on Ubuntu.

`libf2c`: A `fortran` compatibility library that can also be installed on most linux systems using the system package manager

## Mac

`scons`: Can be installed with `homebrew` from the terminal using `brew install scons`

`libf2c`: See instructions [here](http://hpc.sourceforge.net/buildf2c)

# Compilation

`cd` into the folder named **fast** and run `scons` from the command line. The `scons` build system will compile an executable named  `plumeTrackFast` in `build\release\program`.

# Usage

The workflow consists of two parts: (i) CitcomS raw cap files for a given number of timesteps are processed through the `plumeTrackFast` program (ii) the output from the first step are post-processed through a python script to generate a list of plumes detected at each given timestep.

## (i) Running `plumeTrackFast`

`./plumeTrackFast -h` produces the following help message:

```
************************************************************************
* plumeTrackFast (v 0.1)                                               *
************************************************************************

plumeTrackFast features:

USAGE: ./plumeTrackFast [REQUIRED OPTIONS]

OPTIONS:

-h, -help, --help, --usage   Display usage instructions.
--cutoff-percentile ARG      Cut-off percentile below which field values are
                             ignored for cluster analysis.
--data-dir ARG               Path to available cap-files.
--data-file ARG              Model name.
--output-file-basename ARG   Base-name for output files
--radius ARG                 Radius (m)
--start-depth ARG            Depth (m) at which the model-space is scanned for
                             plumes
--start-time ARG             Time (Myr) from which the model-space is scanned
                             for plumes
--stop-time ARG              Time (Myr) from which the model-space is no longer
                             scanned for plumes
--time-file ARG              A two-column text file listing model-times (Myr)
                             and model-time-steps corresponding to available
                             cap-files.
--tracer-flavour ARG         Flavour (1 based) of deep tracers associated with
                             LLSVPs. Set to -1 to ignore tracers.
--validation-depth ARG       Depth (m) at which the model-space is scanned to
                             validate that plume-conduits found at depth
                             '--start-depth' are indeed plumes
--velocity-scale ARG         Velocity-scale (m/Myr)
EXAMPLES:

./plumeTrackFast --data-file <gpm19> --data-dir <pathToCapFiles> --time-file <modelTimes.txt> --velocity-scale 4.953209857165279 --radius 6371e3 --start-time <0> --stop-time <230> --output-file-basename <gpm19.plumes> --start-depth 350e3 --validation-depth 1500e3 --cutoff-percentile 5 --tracer-flavour 5

plumeTrackFast v(0.1), Copyright (C) 2014 Rakib Hassan (rakib.hassan@sydney.edu.au)
This program is free and without warranty.
```

Steps to reproduce the example outputs in `example_output.tar.gz` are as follows:

 * Untar and unzip `data.tar.gz`, which contains cap files for a model at times 0, 65 and 70 Myr. Typically, cap files are output every 5 Myr -- to keep the data volume manageable, data for only three timesteps are provided here. Note that cap files for the 0th timestep are mandetory.
 * create a folder in your work area named e.g. `output`
 * Run the following command to generate cluster analyses output for the above timesteps.
 
 ```
 <path-to-plume-tracker>/fast/build/release/program/plumeTrackFast --data-file gpm58 --data-dir data/gpm58/cap --time-file data/gpm58/cap/times.txt --velocity-scale 4.953209857165279 --radius 6371e3 --start-time 0 --stop-time 230 --output-file-basename output/gpm58.plumes --start-depth 350e3 --validation-depth 1500e3 --cutoff-percentile 5 --tracer-flavour 5
 ```
Note that the `--velocity-scale` parameter is used to convert the non-dimensional model velocities to units of m/Myr and is obtained as described in the CitcomS manual. The above should produce three files, corresponding to the three timesteps, in the output folder. The column-descripts are as follows:

```
#column-name  description

r           : normalized radial distance from centre
theta       : latitude (degrees)
phi         : longitude (degrees)
time        : model time (Myr)
magGradVr   : magnitude of the gradient of radial velocity vector
v_r         : radial velocity (m/Myr)
T           : nondimensional temperature
Eta         : viscosity
data        : magGradVr after applying percentile cutoff
cid         : cluster-id, assigned by clustering algorithm
magGradVrVD : as magGradVr, but at validation-depth
v_rVD       : as v_r, but at validation-depth
TVD         : as T, but at validation-depth
EtaVD       : as Eta, but at validation-depth
dataVD      : as data, but at validation-depth
cidVD       : as cid, but at validation-depth
tracComp    : composition of tracer id given by --tracer-flavour, 
              at a model depth given by --start-depth
```

## (ii) Postprocessing with `plotPlumeAscentRate.py`



