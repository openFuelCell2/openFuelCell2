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

\*---------------------------------------------------------------------------*/

#include "noneDWModel.H"
#include "fuelCellSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dissolvedModels
{
    defineTypeNameAndDebug(noneDWModel, 0);

    addToRunTimeSelectionTable
    (
        dissolvedModel,
        noneDWModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissolvedModels::noneDWModel::noneDWModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dissolvedModel(mesh),

    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dissolvedModels::noneDWModel::~noneDWModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::dissolvedModels::noneDWModel::solve()
{
    // nothing
}


void Foam::dissolvedModels::noneDWModel::correct()
{
    // nothing
}


void Foam::dissolvedModels::noneDWModel::update
(
    const word& clName
)
{
    // nothing
}


void Foam::dissolvedModels::noneDWModel::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    // nothing
}


void Foam::dissolvedModels::noneDWModel::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    // nothing
}


bool Foam::dissolvedModels::noneDWModel::read(const dictionary& dict)
{
    return dict.isDict(type() + "Coeffs");
}
// ************************************************************************* //
