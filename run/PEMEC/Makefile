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
	(cp -rf ./system/controlDict.mesh ./system/controlDict)
	(blockMesh)
	(renumberMesh -overwrite)
	(./preprocessing.sh | tee log.pre)
	(cp -rf ./system/controlDict.run ./system/controlDict)

salomeMesh:
	(cp -rf ./system/controlDict.mesh ./system/controlDict)
	(salome -t salomeMesh.py)
	(ideasUnvToFoam salome.unv)
	(renumberMesh -overwrite)
	(./preprocessing.sh | tee log.pre)
	(cp -rf ./system/controlDict.run ./system/controlDict)

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
	(foamToVTK -latestTime -region air);
	(foamToVTK -latestTime -region fuel);
	(foamToVTK -latestTime -region interconnect);
	(foamToVTK -latestTime -region electrolyte);
	(foamToVTK -latestTime -region phiEC);
	(foamToVTK -latestTime -region phiEA);
	(foamToVTK -latestTime -region phiI);

################################################################################

viewAll:
	(foamToVTK);
	(foamToVTK -region air);
	(foamToVTK -region fuel);
	(foamToVTK -region interconnect);
	(foamToVTK -region electrolyte);
	(foamToVTK -region phiEC);
	(foamToVTK -region phiEA);
	(foamToVTK -region phiI);

################################################################################
