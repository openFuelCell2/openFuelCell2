#!/bin/bash

#----------------------------------------------------------------------#
# Contributor |   Shidong Zhang                                        #
# Email       |   s.zhang@fz-juelich.de                                #
#----------------------------------------------------------------------#
# Solver      |   openFuelCellFoam2                                    #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #
#----------------------------------------------------------------------#
# Source code |   https://jugit.fz-juelich.de/s.zhang/fuelcellfoam     #
# Update from |   15.06.2022                                           #
#----------------------------------------------------------------------#

# Invoke python script
python3 parameterScript.py

cp system/controlDict.mesh system/controlDict

# create cellID
decomposeParID -region air -coordinate 
decomposeParID -region fuel -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'
decomposeParID -region electrolyte -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'
decomposeParID -region interconnect -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'
decomposeParID -region phiEA -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'
decomposeParID -region phiEC -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'
decomposeParID -region phiI -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'
decomposeParID -coordinate '((0.00 -0.001 0) (0.04 0.001 0))'


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

rm processors$NPROCS/constant/polyMesh/cellZones

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.parallel.afei -constant -noZero -parallel -fileHandler collated
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

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.parallel.phiE -constant -noZero -parallel -fileHandler collated
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

mpirun -np $NPROCS topoSet -dict ./system/topoSetDict.parallel.phiI -constant -noZero -parallel -fileHandler collated
mpirun -np $NPROCS splitMeshRegions -cellZonesOnly -parallel -fileHandler collated

# sleep to wait for files
while [! -d processors$NPROCS/1];
do
    sleep 1s
done

cp -rf processors$NPROCS/1/phiI/polyMesh processors$NPROCS/constant/phiI
rm -rf processors$NPROCS/1

rm processors$NPROCS/constant/polyMesh/cellZones

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

cp system/controlDict.run system/controlDict")
