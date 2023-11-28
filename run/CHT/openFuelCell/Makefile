#----------------------------------------------------------------------#
# Solver      |   openFuelCell                                         #
# OpenFOAM    |   OpenFOAM-10                                          #
#----------------------------------------------------------------------#
# Source code |   https://github.com/openFuelCell2/openFuelCell2       #
# Update from |   28.11.2023                                           #
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
	(rm -rf [01]);
	(cp -rf ./system/controlDict.run ./system/controlDict);
	(cp -r 0.orig 0);

################################################################################

decompose:
	(decomposeParID -allRegions -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')

################################################################################

parallel:

	(./parallel.sh | tee log.parallel);

################################################################################

run:
	(mpirun -np 4 openFuelCell -parallel | tee log.run);

################################################################################

plot:
	(gnuplot ResidualPlot)

################################################################################

srun:
	(openFuelCell | tee log.run);

################################################################################

reconstruct:
	(reconstructPar);
	(reconstructPar -allRegions);

################################################################################

view:
	(foamToVTK -latestTime);
	(foamToVTK -latestTime -region hot);
	(foamToVTK -latestTime -region cold);
	(foamToVTK -latestTime -region solid);

################################################################################

viewAll:
	(foamToVTK);
	(foamToVTK -region hot);
	(foamToVTK -region cold);
	(foamToVTK -region solid);

################################################################################
