#!/bin/csh

#----------------------------------------------------------------------#
# Solver      |   chtMultiRegionFoam                                   #
# OpenFOAM    |   OpenFOAM-10                                          #
#----------------------------------------------------------------------#
# Source code |   OpenFOAM-10                                          #
# Update from |   28.11.2023                                           #
#----------------------------------------------------------------------#

## edit system/decomposeParDict for the desired decomposition
## set environment variable NPROCS to number of processors.
##     e.g., setenv NPROCS 2
## make mesh
##
## then:

echo "NPROCS = " $NPROCS

cp system/controlDict.mesh system/controlDict

cp system/decomposeParDict system/hot/.
cp system/decomposeParDict system/cold/.
cp system/decomposeParDict system/solid/.

decomposePar -region hot -fileHandler collated 
decomposePar -region cold -fileHandler collated
decomposePar -region solid -fileHandler collated

cp system/controlDict.run system/controlDict
