/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | openFuelCell
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "diffusivityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusivityModel, 0);
    defineRunTimeSelectionTable(diffusivityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModel::diffusivityModel
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff
)
:
    mesh_(mesh),
    diff_(diff),
    zoneName_(name),
    cells_(),
    firstIndex_(true)
{
    labelList cellZoneIDs = mesh.cellZones().indices(zoneName_);

    if (cellZoneIDs.empty())
    {
        FatalErrorInFunction
            << "cannot find cellZone " << zoneName_
            << exit(FatalError);
    }

    for (auto& id : cellZoneIDs)
    {
        for (auto& cellI : mesh.cellZones()[id])
        {
            cells_.append(cellI);
        }
    }
}

Foam::diffusivityModel::diffusivityModel
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
    const dictionary& dict
)
:
    mesh_(mesh),
    diff_(diff),
    zoneName_(name),
    cells_(),
    dict_(dict),
    firstIndex_(true)
{
    labelList cellZoneIDs = mesh.cellZones().indices(zoneName_);

    if (cellZoneIDs.empty())
    {
        FatalErrorInFunction
            << "cannot find cellZone " << zoneName_
            << exit(FatalError);
    }

    for (auto& id : cellZoneIDs)
    {
        for (auto& cellI : mesh.cellZones()[id])
        {
            cells_.append(cellI);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::diffusivityModel::~diffusivityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// none

// ************************************************************************* //
