/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{
    defineTypeNameAndDebug(linear, 0);
    addToRunTimeSelectionTable
    (
        partitioningModel,
        linear,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
linear::linear(const dictionary& dict)
:
    partitioningModel(),
    alphaLiquid1_(dict.get<scalar>("alphaLiquid1")),
    alphaLiquid0_(dict.get<scalar>("alphaLiquid0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
linear::~linear()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::partitioningModels::
linear::fLiquid
(
    const scalarField& alphaLiquid
) const
{
    return max
        (
            scalar(0),
            min
            (
                scalar(1)-(alphaLiquid1_ - alphaLiquid)
               /(alphaLiquid1_ - alphaLiquid0_),
                scalar(1)
            )
        );
}


void Foam::wallBoilingModels::partitioningModels::
linear::write(Ostream& os) const
{
    partitioningModel::write(os);
    os.writeEntry("alphaLiquid1", alphaLiquid1_);
    os.writeEntry("alphaLiquid0", alphaLiquid0_);
}


// ************************************************************************* //
