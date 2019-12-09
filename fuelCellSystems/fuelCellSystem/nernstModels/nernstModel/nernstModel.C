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

#include "nernstModel.H"
#include "phaseModel.H"

namespace Foam
{
    defineTypeNameAndDebug(nernstModel, 0);
    defineRunTimeSelectionTable(nernstModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::nernstModel::nernstModel
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dictionary& dict
)
:
    volScalarField
    (
        IOobject
        (
            "nernst",
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar
        (
            "nernst",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),

    dict_(dict),
    rxnList_(dict_.lookup("RxnList")),
    deltaH_(this->size(), 0.0),
    deltaS_(deltaH_.size(), 0.0),
    residualY_(dict_.lookupOrDefault<scalar>("residualY", 1.0e-6)),
    pRef_("pRef", dimPressure, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nernstModel::~nernstModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
