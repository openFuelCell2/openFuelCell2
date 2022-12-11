/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    openFuelCell

Description
    Solver for electrochemical devices...

Contributors
    Shidong Zhang (s.zhang@fz-juelich.de)
    Steven B. Beale (s.beale@fz-juelich.de)
    Steffen Hess (s.hess@fz-juelich.de)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "regionTypeList.H"
#include "fuelCellSystem.H"
#include "phaseSystem.H"
#include "regionCourantNo.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    if (!LTS)
    {
        #include "multiRegionCourantNo.H"
        #include "setInitialMultiRegionDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    fuelCell.mapToCell();

    while (runTime.run())
    {
        runTime++;
        Info<< "\nTime = " << runTime.timeName() << nl << endl;

        if (!LTS)
        {
            #include "createTimeControls.H"
            #include "multiRegionCourantNo.H"
            #include "setMultiRegionDeltaT.H"
        }
        else
        {
            fuelCell.setRDeltaT();
        }

        fuelCell.mapFromCell();

        fuelCell.correct();

        fuelCell.solve();

        fuelCell.mapToCell();

        #include "EEqns.H"

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
