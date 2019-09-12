# plume-tracker

Workflow for detecting plumes in CitcomS global mantle-flow models.

# Prerequisites

## Linux
`scons`: A python build utility can be installed on linux systems using the package manager, e.g. `apt` on Ubuntu.

`libf2c`: A `fortran` compatibility library that can also be installed on most linux systems using the system package manager

## Mac

`scons`: Can be installed with `homebrew` as `brew install scons`

`libf2c`: See instructions [here](http://hpc.sourceforge.net/buildf2c)

# Compilation

`cd` into the folder named **fast** and run `scons` from the command line. The `scons` build system will compile an executable named  `plumeTrackFast` in `build\release\program`.

# Usage

The workflow has two parts: (i) CitcomS raw cap files for a given number of time steps are processed through the `plumeTrackFast` program (ii) the output from the first step are post-processed through a python script to generate a list of plumes detected at each given time step.

