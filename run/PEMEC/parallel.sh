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

decomposePar >& log.decompose

cp system/decomposeParDict system/air/.
cp system/decomposeParDict system/fuel/.
cp system/decomposeParDict system/electrolyte/.
cp system/decomposeParDict system/interconnect/.
cp system/decomposeParDict system/phiEA/.
cp system/decomposeParDict system/phiEC/.
cp system/decomposeParDict system/phiI/.

decomposePar -region air
decomposePar -region fuel
decomposePar -region electrolyte
decomposePar -region interconnect
decomposePar -region phiEA
decomposePar -region phiEC
decomposePar -region phiI

# Step 1:
# air/fuel/electrolyte/interconnect

ii=0
while [ -d processor$ii ];
do
  mv processor$ii/constant/polyMesh/cellZones processor$ii/constant/polyMesh/cellZones_bk

  ii=$((ii+1))
done

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.afei -parallel
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel

# sleep to wait for files
while [ ! -d processor0/1 ];
do
    sleep 1s
done

ii=0
while [ -d processor$ii ];
do
  cp -rf processor$ii/1/air/polyMesh processor$ii/constant/air
  cp -rf processor$ii/1/fuel/polyMesh processor$ii/constant/fuel
  cp -rf processor$ii/1/electrolyte/polyMesh processor$ii/constant/electrolyte
  cp -rf processor$ii/1/interconnect/polyMesh processor$ii/constant/interconnect
  rm -rf processor$ii/1

  ii=$((ii+1))
done


# Step 2:
# phiEA, phiEC

rm processor*/constant/polyMesh/cellZones

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.phiE -constant -noZero -parallel
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel

# sleep to wait for files
while [ ! -d processor0/1 ];
do
    sleep 1s
done

ii=0
while [ -d processor$ii ];
do
  cp -rf processor$ii/1/phiEA/polyMesh processor$ii/constant/phiEA
  cp -rf processor$ii/1/phiEC/polyMesh processor$ii/constant/phiEC
  rm -rf processor$ii/1

  ii=$((ii+1))
done

# Step 3:
# phiI

rm processor*/constant/polyMesh/cellZones

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.phiI -constant -noZero -parallel
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel

# sleep to wait for files
while [ ! -d processor0/1 ];
do
    sleep 1s
done

ii=0
while [ -d processor$ii ];
do
  cp -rf processor$ii/1/phiI/polyMesh processor$ii/constant/phiI
  rm -rf processor$ii/1

  mv processor$ii/constant/polyMesh/cellZones_bk processor$ii/constant/polyMesh/cellZones

  ii=$((ii+1))
done


## patches:

# air zones
mpirun -np $NPROCS topoSet -region air -noZero -constant -parallel

# fuel zones
mpirun -np $NPROCS topoSet -region fuel -noZero -constant -parallel

# electric zones
mpirun -np $NPROCS topoSet -region phiEA -noZero -constant -parallel
mpirun -np $NPROCS topoSet -region phiEC -noZero -constant -parallel
mpirun -np $NPROCS topoSet -region phiI -noZero -constant -parallel

rm -rf system/phi0
rm -rf system/phiI0
rm -rf system/phiE0
rm -rf system/phiE1

cp system/controlDict.run system/controlDict
