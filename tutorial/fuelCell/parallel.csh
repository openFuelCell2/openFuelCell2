#!/bin/csh

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

# mpirun -np $NPROCS setSet -batch ./config/make.zoneSet -constant -noVTK -parallel -noZero
mv constant/polyMesh/cellZones constant/polyMesh/cellZones_backup


@ I = 0
while ( $I < $NPROCS )

rm -rf processor$I/constant/*/polyMesh/sets/*

mkdir processor$I/constant/polyMesh/tmp

mv processor$I/constant/polyMesh/sets/* processor$I/constant/polyMesh/tmp/.

cp processor$I/constant/polyMesh/tmp/air processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/fuel processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/electrolyte processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/interconnect processor$I/constant/polyMesh/sets/.

rm processor$I/constant/polyMesh/cellZones

@ I++

end

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -constant -parallel

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec splitMeshRegions -cellZonesOnly -parallel >& log.split1

@ I = 0
while ( $I < $NPROCS )

cp -rf processor$I/1/. processor$I/constant/.
rm -rf processor$I/1

rm processor$I/constant/polyMesh/sets/*

cp processor$I/constant/polyMesh/tmp/phiEA processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/phiEC processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/phiI0 processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/phi0 processor$I/constant/polyMesh/sets/.

rm processor$I/constant/polyMesh/cellZones

@ I++

end

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -constant -parallel

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec splitMeshRegions -cellZonesOnly -parallel >& log.split2

@ I = 0
while ( $I < $NPROCS )

cp -rf processor$I/1/. processor$I/constant/.
rm -rf processor$I/1

rm -rf processor$I/constant/phiI0
rm -rf processor$I/constant/phi0

##mkdir processor$I/constant/polyMesh/sets
rm processor$I/constant/polyMesh/sets/*

cp processor$I/constant/polyMesh/tmp/phiE0 processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/phiE1 processor$I/constant/polyMesh/sets/.
cp processor$I/constant/polyMesh/tmp/phiI processor$I/constant/polyMesh/sets/.

rm processor$I/constant/polyMesh/cellZones

@ I++

end

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -constant -parallel

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec splitMeshRegions -cellZonesOnly -parallel >& log.split3


@ I = 0
while ( $I < $NPROCS )

cp -rf processor$I/1/. processor$I/constant/.

rm -rf processor$I/1

rm -rf processor$I/constant/phiE0
rm -rf processor$I/constant/phiE1

rm processor$I/constant/polyMesh/sets/*

cp processor$I/constant/polyMesh/tmp/* processor$I/constant/polyMesh/sets/.

rm processor$I/constant/polyMesh/cellZones

rm -rf processor$I/constant/polyMesh/tmp

@ I++

end

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -constant -parallel

## patches:

# air zones
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setSet -batch ./config/make.setAir -region air -constant -noVTK -parallel -noZero
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -region air -constant -parallel

# fuel zones
#
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setSet -batch ./config/make.setFuel -region fuel -constant -noVTK -parallel -noZero
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -region fuel -constant -parallel

# electronic zones
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setSet -batch ./config/make.setPhiEA -region phiEA -constant -noVTK -parallel -noZero
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -region phiEA -constant -parallel

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setSet -batch ./config/make.setPhiEC -region phiEC -constant -noVTK -parallel -noZero
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -region phiEC -constant -parallel

$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setSet -batch ./config/make.setPhiI -region phiI -constant -noVTK -parallel -noZero
$MPIEXEC $FLAGS_MPI_BATCH -np $NPROCS foamExec setsToZones -noFlipMap -region phiI -constant -parallel

rm -rf system/phi0
rm -rf system/phiI0
rm -rf system/phiE0
rm -rf system/phiE1

rm -rf VTK
mv constant/polyMesh/cellZones_backup constant/polyMesh/cellZones
cp system/controlDict.run system/controlDict
