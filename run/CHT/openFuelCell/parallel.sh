#!/bin/bash

#----------------------------------------------------------------------#
# Solver      |   openFuelCell                                         #
# OpenFOAM    |   OpenFOAM-10                                          #
#----------------------------------------------------------------------#
# Source code |   https://github.com/openFuelCell2/openFuelCell2       #
# Update from |   28.11.2023                                           #
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

decomposePar

decomposePar -region hot
decomposePar -region cold
decomposePar -region solid

mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel

# sleep to wait for files
ii=0
while [ ! -d processor$ii/1 ];
do
  sleep 1s
  ii=$((ii+1))
done

ii=0
while [ -d processor$ii ];
do
  cp -rf processor$ii/1/hot/polyMesh processor$ii/constant/hot
  cp -rf processor$ii/1/cold/polyMesh processor$ii/constant/cold
  cp -rf processor$ii/1/solid/polyMesh processor$ii/constant/solid

  rm -rf processor$ii/1

  ii=$((ii+1))
done

cp system/controlDict.run system/controlDict
