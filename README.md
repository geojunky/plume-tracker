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

Steps to generate cluster analyses outputs are as follows:

 * In your work area untar and unzip `data.tar.gz`, which contains cap files for a model at times 0, 65 and 70 Myr. Typically, cap files are output every 5 Myr -- to keep the data volume manageable, example data for only three timesteps are provided here. Note that cap files for the 0th timestep are mandatory.
 * create a folder in your work area named e.g. `output`
 * Run the following command to generate cluster analyses output for the above timesteps.
 
 ```
 <path-to-plume-tracker>/fast/build/release/program/plumeTrackFast --data-file gpm58 --data-dir data/gpm58/cap --time-file data/gpm58/cap/times.txt --velocity-scale 4.953209857165279 --radius 6371e3 --start-time 0 --stop-time 230 --output-file-basename output/gpm58.plumes --start-depth 350e3 --validation-depth 1500e3 --cutoff-percentile 5 --tracer-flavour 5
 ```
Note that the `--velocity-scale` parameter is used to convert the non-dimensional model velocities to units of m/Myr and is obtained as described in the CitcomS manual. The above should produce three files, corresponding to the three timesteps, in the output folder. The column-descriptions are as follows:

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
tracComp    : composition (as a fraction) of tracer id  
              given by --tracer-flavour, at a model depth 
              given by --start-depth
```

## (ii) Postprocessing with `plotPlumeAscentRate.py`

`./plotPlumeAscentRate.py -h` produces the following help message:

```
Usage:

Plots results from cluster-analysis from plumeTrackFast

Usage 1: plotPlumeAscentRate.py -b <base file-name> -d <reconstruction-directory> -n <neighbour file-name> -r <radial-velocity minimum> -s <start-age> -o <output file-name>



Options:
  -h, --help            show this help message and exit
  -b BASENAME, --base-name=BASENAME
                        File name for clustered data
  -d RECONSPATH, --directory=RECONSPATH
                        Path to reconstruction files
  -s STARTAGE, --start-age=STARTAGE
                        Start age
  -n NBFILENAME, --neighbours=NBFILENAME
                        File containing node-neighbours from spherical
                        triangulation
  -r RADIALVELOCITYMIN, --radial-velocity-minimum (cm/yr)=RADIALVELOCITYMIN
                        Minimum radial velocity that a plume-conduit must meet
                        to be considered as such
  -o OUTFILENAME, --output-file=OUTFILENAME
                        Output-file name
```

Steps for postprocessing output from (i) are as follows:

* Create another sibling output folder, e.g. `plume350`
* Run the following:

```
<path-to-plume-tracker>/postprocess/plotPlumeAscentRate.py -b ../output/gpm58.plumes.clustered -d ../data/reconstructedShapeFiles/ -n ../data/neighbours.txt -r 10 -s 230 -o plume350/plumes.gpm58.txt
```
Note that the `reconstructedShapeFiles` folder contains reconstructed shapefiles for ages corresponding to the available model times in this example. The `neighbours.txt` contains natural neighbour indices for each node of the CitcomS mesh at the surface and is obtained through a spherical triangulation. Note that this file corresponds to a CitcomS mesh of resolution (129x129x65 x 12); a different resolution mesh would require regenerating this file.

The above should produce the following in the `plume350` folder:

 1. `plumes.gpm58.txt`, a text file (column descriptions below) containing entries, grouped by model time, for plumes detected
 ```
           #column-name        description
          time            :   model time (Myr)
          theta           :   latitude (degrees)
          phi             :   longitude (degrees)
          T               :   nondimensional temperature
          Vr              :   radial velocity (m/Myr)
          Flux            :   buoyancy flux (Mg/s)
          Area            :   conduit area (km2)
          Eruption        :   if a plume wasn't detected within 500 km of 
                              the current location in the preceding 5 Myr, 
                              a detected plume is considered a new eruption
                              and this flag is set to 1; 0 otherwise
          meanBgTemp      :   a secondary mean nondimensional temperature, 
                              computed based on nodes that are warmer than
                              the mean temperature within a ~400 km region
                              around a plume conduit
          avgConduitTemp  :   mean temperature within a plume conduit. 
                              See appendix A in Hassan et al. (2015) for a 
                              detailed account of how plume conduits are
                              delineated
          concentration   :   concentration of anomalous material, as a 
                              fraction, near the surface
 ```
 2. Three plots, corresponding to the model timesteps, showing plumes detected in each
 3. `plume350/supp/`, a folder containing supplementary plots for each detected plume, shown in a local coordinate system
 
