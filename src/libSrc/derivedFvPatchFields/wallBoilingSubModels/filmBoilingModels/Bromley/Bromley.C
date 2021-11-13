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

#include "Bromley.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{
    defineTypeNameAndDebug(Bromley, 0);
    addToRunTimeSelectionTable
    (
        filmBoilingModel,
        Bromley,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::filmBoilingModels::Bromley::Bromley
(
    const dictionary& dict
)
:
    filmBoilingModel(),
    Cn_(dict.getOrDefault<scalar>("Cn", 0.62)),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 1)),
    L_(dict.get<scalar>("L"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::filmBoilingModels::Bromley::~Bromley()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::filmBoilingModels::Bromley::htcFilmBoil
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{

    const fvPatchScalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];
    const uniformDimensionedVectorField& g =
        liquid.mesh().time().lookupObject<uniformDimensionedVectorField>("g");

    const fvPatchScalarField& rhoVaporw
    (
        vapor.thermo().rho()().boundaryField()[patchi]
    );

    const scalarField rhoLiq(liquid.thermo().rho(patchi));
    const scalarField kappaVapor(vapor.kappa(patchi));

    tmp<volScalarField> tCp = vapor.thermo().Cp();
    const volScalarField& Cp = tCp();
    const scalarField& CpVapor = Cp.boundaryField()[patchi];

    const scalarField muVapor(vapor.mu(patchi));
    //const scalarField dbVapor(vapor.d()().boundaryField()[patchi]);

    const scalarField htcRad
    (
        emissivity_*physicoChemical::sigma.value()*(pow4(Tw) - pow4(Tsatw))
      / max((Tw - Tsatw), scalar(1e-4))
    );

    return
        Cn_*pow
        (
            pow3(kappaVapor)
           *rhoVaporw*(rhoLiq - rhoVaporw)*mag(g.value())
           *(L + 0.4*CpVapor*max((Tw-Tsatw), scalar(0)))
           /(L_*muVapor*max((Tw-Tsatw), scalar(1e-4))),
            0.25
        ) + 0.75*htcRad;
}


void Foam::wallBoilingModels::filmBoilingModels::Bromley::write
(
    Ostream& os
) const
{
    filmBoilingModel::write(os);
    os.writeEntry("Cn", Cn_);
    os.writeEntry("L", L_);
    os.writeEntry("emissivity", emissivity_);
}


// ************************************************************************* //
