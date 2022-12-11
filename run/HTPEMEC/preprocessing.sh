/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

# Rename the initial field
# mv 0 0.orig

SECONDS=0

# Step 1:
# air/fuel/electrolyte/interconnect

topoSet -dict ./system/topoSetDict.serial.afei -noZero -constant
splitMeshRegions -cellZonesOnly

# copy
cp -r 1/air/polyMesh constant/air/.
cp -r 1/fuel/polyMesh constant/fuel/.
cp -r 1/electrolyte/polyMesh constant/electrolyte/.
cp -r 1/interconnect/polyMesh constant/interconnect/.

rm -rf 1

topoSet -region air
topoSet -region fuel

# Step 2:
# phiEA/phiEC

rm constant/polyMesh/cellZones
topoSet -dict ./system/topoSetDict.serial.phiE -noZero -constant
splitMeshRegions -cellZonesOnly

cp -r 1/phiEC/polyMesh constant/phiEC/.
cp -r 1/phiEA/polyMesh constant/phiEA/.

rm -rf 1

topoSet -region phiEC
topoSet -region phiEA

rm -rf constant/phiI0
rm -rf constant/phi0
rm -rf system/phiI0
rm -rf system/phi0
rm -rf 0/phiI0
rm -rf 0/phi0

# Step 3:
# phiI

rm constant/polyMesh/cellZones
topoSet -dict ./system/topoSetDict.serial.phiI -noZero -constant
splitMeshRegions -cellZonesOnly

cp -r 1/phiI/polyMesh constant/phiI/.

rm -rf 1

topoSet -region phiI

rm -rf constant/phiE1
rm -rf constant/phiE0
rm -rf system/phiE1
rm -rf system/phiE0
rm -rf 0/phiE1
rm -rf 0/phiE0

rm constant/polyMesh/cellZones
topoSet -dict ./system/topoSetDict.parallel.afei -noZero -constant

# Rename the original field to 0
#rm -rf 0
#mv 0.orig 0

duration=$SECONDS
echo -e "  -Executation time for preprocessing is : $duration Sec \n\n"
