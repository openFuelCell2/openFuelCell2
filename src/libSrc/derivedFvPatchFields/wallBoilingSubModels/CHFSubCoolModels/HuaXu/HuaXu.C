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

#include "HuaXu.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phasePairKey.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{
    defineTypeNameAndDebug(HuaXu, 0);
    addToRunTimeSelectionTable
    (
        CHFSubCoolModel,
        HuaXu,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::CHFModels::HuaXu::HuaXu
(
    const dictionary& dict
)
:
    CHFSubCoolModel(),
    Kburn_(dict.getOrDefault<scalar>("Kburn", 1.5))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::CHFModels::HuaXu::~HuaXu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::CHFModels::HuaXu::CHFSubCool
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    const uniformDimensionedVectorField& g =
        liquid.mesh().time().lookupObject<uniformDimensionedVectorField>("g");

    const scalarField alphaLiq(liquid.alpha(patchi));

    const scalarField rhoVapor(vapor.thermo().rho(patchi));
    const scalarField rhoLiq(liquid.thermo().rho(patchi));
    tmp<volScalarField> tCp = liquid.thermo().Cp();
    const volScalarField& Cp = tCp();
    const fvPatchScalarField& Cpw = Cp.boundaryField()[patchi];

    const phasePairKey pair(liquid.name(), vapor.name());
    const scalarField sigma
    (
        liquid.fluid().sigma(pair)().boundaryField()[patchi]
    );

    const scalarField Pe
    (
        pow(sigma, 0.75)
       /
        (
            alphaLiq
          * pow(mag(g.value())*(rhoLiq-rhoVapor), 0.25)
          * sqrt(rhoVapor)
        )
    );

    const scalarField Ja
    (
        rhoLiq*Cpw*max(Tsatw - Tl, scalar(0))/(rhoVapor*L)
    );

    return
        Kburn_*(1 + 0.345*Ja/pow(Pe, 0.25));
}


void Foam::wallBoilingModels::CHFModels::HuaXu::write
(
    Ostream& os
) const
{
    CHFSubCoolModel::write(os);
    os.writeEntry("Kburn", Kburn_);
}


// ************************************************************************* //
