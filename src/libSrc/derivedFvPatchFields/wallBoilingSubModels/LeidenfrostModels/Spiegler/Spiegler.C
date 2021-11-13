/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "Spiegler.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace LeidenfrostModels
{
    defineTypeNameAndDebug(Spiegler, 0);
    addToRunTimeSelectionTable
    (
        LeidenfrostModel,
        Spiegler,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::LeidenfrostModels::Spiegler::Spiegler
(
    const dictionary& dict
)
:
    LeidenfrostModel(),
    Tcrit_(dict.getOrDefault<scalar>("Tcrit", 374))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::LeidenfrostModels::Spiegler::~Spiegler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::LeidenfrostModels::Spiegler::TLeid
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            liquid.thermo().p().boundaryField()[patchi].size(),
            27*Tcrit_/32
        )
    );
}


void Foam::wallBoilingModels::LeidenfrostModels::Spiegler::write
(
    Ostream& os
) const
{
    LeidenfrostModel::write(os);
    os.writeEntry("Tcrit", Tcrit_);
}

// ************************************************************************* //
