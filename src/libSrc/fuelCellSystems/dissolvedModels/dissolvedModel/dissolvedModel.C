    /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "dissolvedModel.H"
#include "zeroGradientFvPatchFields.H"

const Foam::word Foam::dissolvedModel::modelName("dissolved");

namespace Foam
{
    defineTypeNameAndDebug(dissolvedModel, 0);
    defineRunTimeSelectionTable(dissolvedModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissolvedModel::dissolvedModel
(
    const fvMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            modelName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    lambda_
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    dmdt_
    (
        IOobject
        (
            "dmdt",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "dmdt",
            dimensionSet(0, -3, -1, 0, 1, 0, 0),
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    Dwm_
    (
        IOobject
        (
            "Dwm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Dwm", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1.0e-6),
        zeroGradientFvPatchScalarField::typeName
    ),

    act_
    (
        IOobject
        (
            "act",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("act", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dissolvedModel::~dissolvedModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dissolvedModel::writeData(Ostream& os) const
{
    return true;
}
