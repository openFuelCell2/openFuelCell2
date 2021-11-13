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

#include "Jeschar.H"
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
    defineTypeNameAndDebug(Jeschar, 0);
    addToRunTimeSelectionTable
    (
        MHFModel,
        Jeschar,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::CHFModels::Jeschar::Jeschar
(
    const dictionary& dict
)
:
    MHFModel(),
    Kmhf_(dict.getOrDefault<scalar>("Kmhf", 1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::CHFModels::Jeschar::~Jeschar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::CHFModels::Jeschar::MHF
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

    const scalarField rhoVapor(vapor.thermo().rho(patchi));
    const scalarField rhoLiq(liquid.thermo().rho(patchi));

    const phasePairKey pair(liquid.name(), vapor.name());
    const scalarField sigma
    (
        liquid.fluid().sigma(pair)().boundaryField()[patchi]
    );

    return
        Kmhf_*0.09*rhoVapor*L
       *(
            pow(sigma/(mag(g.value())*(rhoLiq - rhoVapor)), 0.25)
          * sqrt(mag(g.value())*(rhoLiq - rhoVapor)/(rhoLiq + rhoVapor))
        );
}


void Foam::wallBoilingModels::CHFModels::Jeschar::write
(
    Ostream& os
) const
{
    MHFModel::write(os);
    os.writeEntry("Kmhf", Kmhf_);
}


// ************************************************************************* //
