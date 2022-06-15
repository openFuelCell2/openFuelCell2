/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | openFuelCell
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
    decomposePar IDs

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "regionProperties.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "coordinate",
        "(pointA(min) pointB(max))",
        "Decompose the mesh(block, cube) between <pointA> and <pointB> "
        "- eg, '( (-1 -1 0) (1 1 0) )'"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary cellProperties
    (
        IOobject
        (
            "cellProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label nx(cellProperties.get<label>("nx"));
    label ny(cellProperties.get<label>("ny"));

    // Get all region names
    wordList regionNames;
    if (args.found("allRegions"))
    {
        regionNames = regionProperties(runTime).names();

        Info<< "Decomposing all regions in regionProperties" << nl
            << "    " << flatOutput(regionNames) << nl << endl;
    }
    else
    {
        regionNames.resize(1);
        regionNames.first() =
            args.getOrDefault<word>("region", fvMesh::defaultRegion);
    }

    forAll(regionNames, regioni)
    {

        Foam::word regionName = regionNames[regioni];

        Foam::fvMesh mesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );

        labelIOList meshCellID
        (
            IOobject
            (
                "cellID",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            labelList(mesh.nCells(), 1)
        );

        const bool doDecomposeBlock = args.found("coordinate");

        // this is not actually stringent enough:
        if (!doDecomposeBlock)
        {
            FatalErrorInFunction
                << "No options supplied, please indicate decompose "
                << "block coordinates. "
                << exit(FatalError);
        }

        Pair<vector> c1c2
        (
            args.lookup("coordinate")()
        );

        label NpX(0), NpY(0);

        scalar lengthX = mag(c1c2[0][0] - c1c2[1][0]);
        scalar lengthY = mag(c1c2[0][1] - c1c2[1][1]);

        scalar smallX = min(c1c2[0][0], c1c2[1][0]);
        scalar smallY = min(c1c2[0][1], c1c2[1][1]);

        forAll(meshCellID, cellI)
        {
            scalar cx = mesh.C()[cellI].component(0);
            scalar cy = mesh.C()[cellI].component(1);

            NpX = label((cx - smallX)/lengthX*nx);
            NpY = label((cy - smallY)/lengthY*ny);

            NpX = min(NpX, nx - 1);
            NpY = min(NpY, ny - 1);

            NpX = max(NpX, 0);
            NpY = max(NpY, 0);

            meshCellID[cellI] = nx*NpY + NpX;
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        Info<< "Writing cellID into directory " << meshCellID.path() << nl << endl;
        meshCellID.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
