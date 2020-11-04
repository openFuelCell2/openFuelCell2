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

#include "Schroeder.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::R;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{
    defineTypeNameAndDebug(Schroeder, 0);
    addToRunTimeSelectionTable
    (
        TDNBModel,
        Schroeder,
        dictionary
    );
}
}
}

using Foam::constant::physicoChemical::R;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::TDNBModels::Schroeder::Schroeder
(
    const dictionary& dict
)
:
    TDNBModel(),
    kg_(dict.getOrDefault<scalar>("kg", 1.666))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::TDNBModels::Schroeder::~Schroeder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::TDNBModels::Schroeder::TDNB
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    // Converting from g/mol to  Kg/mol
    const scalarField W(1e-3*liquid.thermo().W()().boundaryField()[patchi]);

    // isoentropic expansion factor for ideal gases

    return
        Tsatw
        /
        (
            1 - log(2*kg_ + 1)*(R.value()*Tsatw)/(W*L)
        );
}


void Foam::wallBoilingModels::TDNBModels::Schroeder::write
(
    Ostream& os
) const
{
    TDNBModel::write(os);
    os.writeEntry("kg", kg_);
}

// ************************************************************************* //
