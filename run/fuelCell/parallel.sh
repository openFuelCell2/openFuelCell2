#!/bin/bash

## edit system/decomposeParDict for the desired decomposition
## set environment variable NPROCS to number of processors.
##     e.g., setenv NPROCS 2
## make mesh
##
## then:
##!/bin/csh

## edit system/decomposeParDict for the desired decomposition
## set environment variable NPROCS to number of processors.
##     e.g., setenv NPROCS 2
## make mesh
##
## then:

echo "NPROCS = " $NPROCS

# To reconstruct and visualize the regions, we need the *ProcAddressing files
# created by decomposePar -region <region name>
# After the region decomp, we rename the processor* directories as proc_*
# to (a) allow the parallel decomp to proceed 
# while (b) saving the *ProcAddressing files for later copy
#

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

# mpirun -np $NPROCS setSet -batch ./config/make.zoneSet -constant -noVTK -parallel -noZero
mv constant/polyMesh/cellZones constant/polyMesh/cellZones_backup

rm -rf processors$NPROCS/constant/*/polyMesh/sets/*

mkdir processors$NPROCS/constant/polyMesh/tmp

mv processors$NPROCS/constant/polyMesh/sets/* processors$NPROCS/constant/polyMesh/tmp/.

## air, fuel, electrolyte, interconnect

cp processors$NPROCS/constant/polyMesh/tmp/air          processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/fuel         processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/electrolyte  processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/interconnect processors$NPROCS/constant/polyMesh/sets/.

rm processors$NPROCS/constant/polyMesh/cellZones

mpirun -np $NPROCS setsToZones -noFlipMap -constant -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

rm processors$NPROCS/constant/polyMesh/sets/*

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
  sleep 1s
done
ls processors$NPROCS/1

cp -rf processors$NPROCS/1/. processors$NPROCS/constant/.
rm -rf processors$NPROCS/1

## phiEA, phiEC

cp processors$NPROCS/constant/polyMesh/tmp/phiEA processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/phiEC processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/phiI0 processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/phi0  processors$NPROCS/constant/polyMesh/sets/.

rm processors$NPROCS/constant/polyMesh/cellZones

mpirun -np $NPROCS setsToZones -noFlipMap -constant -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated 

rm processors$NPROCS/constant/polyMesh/sets/*

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
  sleep 1s
done
ls processors$NPROCS/1

cp -rf processors$NPROCS/1/phiEA processors$NPROCS/constant/.
cp -rf processors$NPROCS/1/phiEC processors$NPROCS/constant/.

rm -rf processors$NPROCS/1

## phiI

cp processors$NPROCS/constant/polyMesh/tmp/phiE0 processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/phiE1 processors$NPROCS/constant/polyMesh/sets/.
cp processors$NPROCS/constant/polyMesh/tmp/phiI  processors$NPROCS/constant/polyMesh/sets/.

rm processors$NPROCS/constant/polyMesh/cellZones

mpirun -np $NPROCS setsToZones -noFlipMap -constant -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

rm processors$NPROCS/constant/polyMesh/sets/*

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
  sleep 1s
done
ls processors$NPROCS/1

cp -rf processors$NPROCS/1/phiI processors$NPROCS/constant/.
rm -rf processors$NPROCS/1

cp processors$NPROCS/constant/polyMesh/tmp/* processors$NPROCS/constant/polyMesh/sets/.

rm processors$NPROCS/constant/polyMesh/cellZones

rm -rf processors$NPROCS/constant/polyMesh/tmp

mpirun -np $NPROCS setsToZones -noFlipMap -constant -parallel -fileHandler collated

## patches:

# air zones
mpirun -np $NPROCS setSet -batch ./config/make.setAir -region air -constant -noVTK -parallel -noZero -fileHandler collated
mpirun -np $NPROCS setsToZones -noFlipMap -region air -constant -parallel -fileHandler collated

# fuel zones
#
mpirun -np $NPROCS setSet -batch ./config/make.setFuel -region fuel -constant -noVTK -parallel -noZero -fileHandler collated
mpirun -np $NPROCS setsToZones -noFlipMap -region fuel -constant -parallel -fileHandler collated

# electronic zones
mpirun -np $NPROCS setSet -batch ./config/make.setPhiEA -region phiEA -constant -noVTK -parallel -noZero -fileHandler collated
mpirun -np $NPROCS setsToZones -noFlipMap -region phiEA -constant -parallel -fileHandler collated

mpirun -np $NPROCS setSet -batch ./config/make.setPhiEC -region phiEC -constant -noVTK -parallel -noZero -fileHandler collated
mpirun -np $NPROCS setsToZones -noFlipMap -region phiEC -constant -parallel -fileHandler collated

mpirun -np $NPROCS setSet -batch ./config/make.setPhiI -region phiI -constant -noVTK -parallel -noZero -fileHandler collated
mpirun -np $NPROCS setsToZones -noFlipMap -region phiI -constant -parallel -fileHandler collated

rm -rf system/phi0
rm -rf system/phiI0
rm -rf system/phiE0
rm -rf system/phiE1

rm -rf VTK
mv constant/polyMesh/cellZones_backup constant/polyMesh/cellZones
cp system/controlDict.run system/controlDict
