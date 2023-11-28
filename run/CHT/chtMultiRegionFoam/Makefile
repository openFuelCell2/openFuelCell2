#----------------------------------------------------------------------#
# Solver      |   chtMultiRegionSimpleFoam                             #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #
#----------------------------------------------------------------------#
# Source code |   https://jugit.fz-juelich.de/s.zhang/fuelcellfoam     #
# Update from |   14.09.2023                                           #
#----------------------------------------------------------------------#

clear:
	(./Allclean)

################################################################################

mesh:
	(cp -rf ./system/controlDict.mesh ./system/controlDict);
	(blockMesh);
	(splitMeshRegions -cellZonesOnly);
	(cp -rf 1/hot/polyMesh constant/hot/.);
	(cp -rf 1/cold/polyMesh constant/cold/.);
	(cp -rf 1/solid/polyMesh constant/solid/.);
	(rm -rf 1);
	(cp -rf ./system/controlDict.run ./system/controlDict);

################################################################################

decompose:
	(decomposeParID -allRegions -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')

################################################################################

parallel:

	(./parallel.sh | tee log.parallel);

################################################################################

run:
	(mpirun -np 4 chtMultiRegionFoam -parallel | tee log.run);

################################################################################

plot:
	(gnuplot ResidualPlot)

################################################################################

srun:
	(chtMultiRegionFoam | tee log.run);

################################################################################

reconstruct:
	(reconstructPar -region hot);
	(reconstructPar -region solid);
	(reconstructPar -region cold);

################################################################################

view:
	(foamToVTK -latestTime -region hot);
	(foamToVTK -latestTime -region cold);
	(foamToVTK -latestTime -region solid);

################################################################################

viewAll:
	(foamToVTK -region hot);
	(foamToVTK -region cold);
	(foamToVTK -region solid);

################################################################################
