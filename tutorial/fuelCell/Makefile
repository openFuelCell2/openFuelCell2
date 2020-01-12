################################################################################
###############
clear:
	(./Allclean)
###############
mesh:
	(cp -rf ./system/controlDict.mesh ./system/controlDict);
	(blockMesh);
	(./pre | tee log.pre);
	(cp -rf ./system/controlDict.run ./system/controlDict);
#############
parallel:
	(./parallel.csh | tee log.parallel);
#############
run:
	(mpirun -np 3 fuelCell0Foam -parallel | tee log.prun);
pPlot:
	(gnuplot pResidualPlot)
srun:
	(fuelCell0Foam | tee log.srun);
sPlot:
	(gnuplot sResidualPlot)
#############
###reconstructMesh:
############
reconstruct:
	(reconstructParMesh);
	(reconstructParMesh -region air 					-constant);
	(reconstructParMesh -region fuel 					-constant);
	(reconstructParMesh -region electrolyte 	-constant);
	(reconstructParMesh -region interconnect 	-constant);
	(reconstructParMesh -region phiEA 				-constant);
	(reconstructParMesh -region phiEC 				-constant);
	(reconstructParMesh -region phiI				 -constant);
	(reconstructPar);
	(reconstructPar -region air);
	(reconstructPar -region fuel);
	(reconstructPar -region electrolyte);
	(reconstructPar -region interconnect);
	(reconstructPar -region phiEA);
	(reconstructPar -region phiEC);
	(reconstructPar -region phiI);
#############
view:
	(foamToVTK -latestTime);
	(foamToVTK -latestTime -region air);
	(foamToVTK -latestTime -region fuel);
	(foamToVTK -latestTime -region electrolyte);
	(foamToVTK -latestTime -region interconnect);
	(foamToVTK -latestTime -region phiEA);
	(foamToVTK -latestTime -region phiEC);
	(foamToVTK -latestTime -region phiI);
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
