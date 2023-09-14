#!/bin/bash

#----------------------------------------------------------------------#
# Solver      |   openFuelCell                                         #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #
#----------------------------------------------------------------------#
# Source code |   https://github.com/openFuelCell2/openFuelCell2       #
# Update from |   14.09.2023                                           #
#----------------------------------------------------------------------#

# Rename the original field to 0
rm -rf 0
cp -r 0.orig 0

SECONDS=0

topoSet -dict ./system/topoSetDict.zoneToSet -noZero -constant

# Step 1:
# air/fuel/electrolyte/interconnect
# backup the cellZones
mv constant/polyMesh/cellZones constant/polyMesh/cellZones_bk
topoSet -dict ./system/topoSetDict.afei -noZero -constant
splitMeshRegions -cellZonesOnly

# copy
cp -r 1/air/polyMesh constant/air/.
cp -r 1/fuel/polyMesh constant/fuel/.
cp -r 1/electrolyte/polyMesh constant/electrolyte/.
cp -r 1/interconnect/polyMesh constant/interconnect/.

rm -rf 1

# Step 2:
# phiEA/phiEC

rm constant/polyMesh/cellZones
topoSet -dict ./system/topoSetDict.phiE -noZero -constant
splitMeshRegions -cellZonesOnly

cp -r 1/phiEC/polyMesh constant/phiEC/.
cp -r 1/phiEA/polyMesh constant/phiEA/.

rm -rf 1
rm -rf constant/phiE0
rm -rf system/phiE0
rm -rf 0/phiE0

# Step 3:
# phiI

rm constant/polyMesh/cellZones
topoSet -dict ./system/topoSetDict.phiI -noZero -constant
splitMeshRegions -cellZonesOnly

cp -r 1/phiI/polyMesh constant/phiI/.

rm -rf 1
rm -rf constant/phiE0
rm -rf system/phiE0
rm -rf 0/phiE0

# mv back the original cell zones
mv constant/polyMesh/cellZones_bk constant/polyMesh/cellZones

topoSet -region air
topoSet -region fuel
topoSet -region phiEC
topoSet -region phiEA
topoSet -region phiI

# Rename the original field to 0
#rm -rf 0
#mv 0.orig 0

duration=$SECONDS
echo -e "  -Executation time for preprocessing is : $duration Sec \n\n"
