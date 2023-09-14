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

decomposePar -fileHandler collated >& log.decompose

cp system/decomposeParDict system/air/.
cp system/decomposeParDict system/fuel/.
cp system/decomposeParDict system/electrolyte/.
cp system/decomposeParDict system/interconnect/.
cp system/decomposeParDict system/phiEA/.
cp system/decomposeParDict system/phiEC/.
cp system/decomposeParDict system/phiI/.

decomposePar -region air -fileHandler collated 
decomposePar -region fuel -fileHandler collated
decomposePar -region electrolyte -fileHandler collated
decomposePar -region interconnect -fileHandler collated
decomposePar -region phiEA -fileHandler collated
decomposePar -region phiEC -fileHandler collated
decomposePar -region phiI -fileHandler collated

# Step 1:
# air/fuel/electrolyte/interconnect

mv processors$NPROCS/constant/polyMesh/cellZones processors$NPROCS/constant/polyMesh/cellZones_bk

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.afei -constant -noZero -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
    sleep 1s
done

cp -rf processors$NPROCS/1/air/polyMesh processors$NPROCS/constant/air
cp -rf processors$NPROCS/1/fuel/polyMesh processors$NPROCS/constant/fuel
cp -rf processors$NPROCS/1/electrolyte/polyMesh processors$NPROCS/constant/electrolyte
cp -rf processors$NPROCS/1/interconnect/polyMesh processors$NPROCS/constant/interconnect

rm -rf processors$NPROCS/1

# Step 2:
# phiEA, phiEC

rm processors$NPROCS/constant/polyMesh/cellZones

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.phiE -constant -noZero -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
    sleep 1s
done

cp -rf processors$NPROCS/1/phiEA/polyMesh processors$NPROCS/constant/phiEA
cp -rf processors$NPROCS/1/phiEC/polyMesh processors$NPROCS/constant/phiEC

rm -rf processors$NPROCS/1

# Step 3:
# phiI

rm processors$NPROCS/constant/polyMesh/cellZones

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.phiI -constant -noZero -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
    sleep 1s
done

cp -rf processors$NPROCS/1/phiI/polyMesh processors$NPROCS/constant/phiI
rm -rf processors$NPROCS/1

mv processors$NPROCS/constant/polyMesh/cellZones_bk processors$NPROCS/constant/polyMesh/cellZones

## patches:

# air zones
mpirun -np $NPROCS topoSet -region air -noZero -constant -parallel -fileHandler collated

# fuel zones
mpirun -np $NPROCS topoSet -region fuel -noZero -constant -parallel -fileHandler collated

# electric zones
mpirun -np $NPROCS topoSet -region phiEA -noZero -constant -parallel -fileHandler collated
mpirun -np $NPROCS topoSet -region phiEC -noZero -constant -parallel -fileHandler collated
mpirun -np $NPROCS topoSet -region phiI -noZero -constant -parallel -fileHandler collated

rm -rf system/phi0
rm -rf system/phiI0
rm -rf system/phiE0
rm -rf system/phiE1

cp system/controlDict.run system/controlDict
