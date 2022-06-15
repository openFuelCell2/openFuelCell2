################################################################################
                                                                               #
#----------------------------------------------------------------------#       #
# Contributor |                                                        #       #
# Email       |   Shidong Zhang (s.zhang@fz-juelich.de)                #       #
#----------------------------------------------------------------------#       #
# Solver      |   openFuelCell                                         #       #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #       #
#----------------------------------------------------------------------#       #
# Source code |   https://jugit.fz-juelich.de/s.zhang/fuelcellfoam     #       #
# Update from |   15.06.2022                                           #       #
#----------------------------------------------------------------------#       #
                                                                               #
################################################################################

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
	(decomposeParID -region hot -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region cold -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region solid -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')

################################################################################

parallel:

	(./parallel.sh | tee log.parallel);

################################################################################

run:
	(mpirun -np 3 openFuelCell -parallel -fileHandler collated | tee log.run);

################################################################################

plot:
	(gnuplot ResidualPlot)

################################################################################

srun:
	(openFuelCell | tee log.run);

################################################################################

reconstruct:
	(reconstructParMesh);
	(reconstructParMesh -region hot -constant -fileHandler collated);
	(reconstructParMesh -region solid -constant -fileHandler collated);
	(reconstructParMesh -region cold -constant -fileHandler collated);
	(reconstructPar -fileHandler collated);
	(reconstructPar -region hot -fileHandler collated);
	(reconstructPar -region solid -fileHandler collated);
	(reconstructPar -region cold -fileHandler collated);

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
