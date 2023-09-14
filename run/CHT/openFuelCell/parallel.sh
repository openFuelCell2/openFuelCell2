#!/bin/bash

#----------------------------------------------------------------------#
# Solver      |   openFuelCell                                         #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #
#----------------------------------------------------------------------#
# Source code |   https://github.com/openFuelCell2/openFuelCell2       #
# Update from |   14.09.2023                                           #
#----------------------------------------------------------------------#

## edit system/decomposeParDict for the desired decomposition
## set environment variable NPROCS to number of processors.
##     e.g., export NPROCS=2
## make mesh
##
## then:

echo "NPROCS = " $NPROCS

cp system/controlDict.mesh system/controlDict

cp system/decomposeParDict system/hot/.
cp system/decomposeParDict system/cold/.
cp system/decomposeParDict system/solid/.

decomposePar -fileHandler collated

decomposePar -region hot -fileHandler collated
decomposePar -region cold -fileHandler collated
decomposePar -region solid -fileHandler collated

mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

# sleep to wait for files
while [ ! -d processors$NPROCS/1 ]
do
    sleep 1s
done

cp -rf processors$NPROCS/1/hot/polyMesh processors$NPROCS/constant/hot
cp -rf processors$NPROCS/1/cold/polyMesh processors$NPROCS/constant/cold
cp -rf processors$NPROCS/1/solid/polyMesh processors$NPROCS/constant/solid

rm -rf processors$NPROCS/1

cp system/controlDict.run system/controlDict
