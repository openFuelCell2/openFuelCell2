// create fuel cell system with mesh...

Info << "Creating fuel cell system\n" << endl;

fuelCellSystem fuelCell(mesh);

regionTypeList& regions = fuelCell.regions();

// temperature
volScalarField& TCell = fuelCell.T();
volScalarField& kCell = fuelCell.k();
volScalarField& contErrCpCell = fuelCell.contErrCp();
volScalarField& rhoCpCell = fuelCell.rhoCp();
surfaceScalarField& rhoCpPhiCell = fuelCell.rhoCpPhi();
fvScalarMatrix& QdotCell = fuelCell.Qdot();

bool LTS = false;

forAll(regions, regionI)
{
    LTS = LTS || fv::localEulerDdt::enabled(regions[regionI]);
}
