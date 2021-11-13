/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationOverpotentialModel.H"
#include "ActivationOverpotentialModel.H"
#include "ButlerVolmer.H"
#include "Tafel.H"
#include "Constant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReactionThermo.H"

#include "thermoPhysicsTypes.H"

#include "rhoConst.H"
#include "perfectFluid.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeActivationOverpotentialModel(Model, Thermo)           \
                                                                               \
    /* Composition at an interface with a multi-component mixture */           \
    makeActivationOverpotentialType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo                                        \
    );                                                                         \
    makeActivationOverpotentialType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo                                              \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace activationOverpotentialModels;

    makeActivationOverpotentialModel
    (
        ButlerVolmer,
        gasEThermoPhysics
    );
    makeActivationOverpotentialModel
    (
        ButlerVolmer,
        constGasEThermoPhysics
    );

    makeActivationOverpotentialModel
    (
        Tafel,
        gasEThermoPhysics
    );
    makeActivationOverpotentialModel
    (
        Tafel,
        constGasEThermoPhysics
    );

    makeActivationOverpotentialModel
    (
        Constant,
        gasEThermoPhysics
    );
    makeActivationOverpotentialModel
    (
        Constant,
        constGasEThermoPhysics
    );
}

// ************************************************************************* //
