################################################################################
                                                                               #
#----------------------------------------------------------------------#       #
# Contributor |                                                        #       #
# Email       |   Shidong Zhang (s.zhang@fz-juelich.de)                #       #
# Email       |   Steffen Hess  (s.hess@fz-juelich.de)                 #       #
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
	(cp -rf ./system/controlDict.mesh ./system/controlDict)
	(blockMesh)
	(renumberMesh -overwrite)
	(./preprocessing.sh | tee log.pre)
	(cp -rf ./system/controlDict.run ./system/controlDict)

################################################################################

decompose:
	(decomposeParID -region air -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region fuel -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region electrolyte -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region interconnect -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region phiEA -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region phiEC -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -region phiI -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')
	(decomposeParID -coordinate '((0.00 -0.001 0) (0.04 0.001 0))')

################################################################################

parallel:
	(./parallel.sh | tee log.parallel);

################################################################################

run:
	(mpirun -np 4 openFuelCell -parallel -fileHandler collated | tee log.run);

################################################################################

plot:
	(gnuplot ResidualPlot)

################################################################################

srun:
	(openFuelCell | tee log.run);

################################################################################

reconstruct:
	(reconstructParMesh);
	(reconstructParMesh -region air 					-constant -fileHandler collated);
	(reconstructParMesh -region fuel 					-constant -fileHandler collated);
	(reconstructParMesh -region electrolyte 	-constant -fileHandler collated);
	(reconstructParMesh -region interconnect 	-constant -fileHandler collated);
	(reconstructParMesh -region phiEA 				-constant -fileHandler collated);
	(reconstructParMesh -region phiEC 				-constant -fileHandler collated);
	(reconstructParMesh -region phiI				 	-constant -fileHandler collated);
	(reconstructPar -fileHandler collated);
	(reconstructPar -region air -fileHandler collated);
	(reconstructPar -region fuel -fileHandler collated);
	(reconstructPar -region electrolyte -fileHandler collated);
	(reconstructPar -region interconnect -fileHandler collated);
	(reconstructPar -region phiEA -fileHandler collated);
	(reconstructPar -region phiEC -fileHandler collated);
	(reconstructPar -region phiI -fileHandler collated);

################################################################################

view:
	(foamToVTK -latestTime);
	(foamToVTK -latestTime -region air);
	(foamToVTK -latestTime -region fuel);
	(foamToVTK -latestTime -region electrolyte);
	(foamToVTK -latestTime -region interconnect);
	(foamToVTK -latestTime -region phiEA);
	(foamToVTK -latestTime -region phiEC);
	(foamToVTK -latestTime -region phiI);

################################################################################

viewAll:
	(foamToVTK);
	(foamToVTK -region air);
	(foamToVTK -region fuel);
	(foamToVTK -region electrolyte);
	(foamToVTK -region interconnect);
	(foamToVTK -region phiEA);
	(foamToVTK -region phiEC);
	(foamToVTK -region phiI);

################################################################################
