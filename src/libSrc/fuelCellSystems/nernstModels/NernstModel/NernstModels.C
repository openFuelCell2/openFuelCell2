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

#include "nernstModel.H"
#include "NernstModel.H"
#include "standard.H"
#include "fixedValue.H"

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

#define makeSpecieNernstModel(Model, Thermo1, Thermo2)           \
                                                                               \
    /* Composition at an interface with a multi-component mixture */           \
    makeSpecieNernstType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo2                                         \
    );                                                                         \
    makeSpecieNernstType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo1,                                              \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo2                                         \
    );                                                                         \
                                                                               \
    /* Composition at an interface with a reacting mixture */                  \
    makeSpecieNernstType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo2                                               \
    );                                                                         \
    makeSpecieNernstType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo1,                                              \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo2                                               \
    );

#define makeNernstModel(Model, Thermo1, Thermo2)                 \
                                                                               \
    /* Composition at an interface with a pure mixture */                      \
    makeNernstType                                               \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        multiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoThermo,                                                \
        pureMixture, Thermo2                                                   \
    );                                                                         \
    makeNernstType                                               \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        reactingMixture, Thermo1,                                              \
        heRhoThermo, rhoThermo,                                                \
        pureMixture, Thermo2                                                   \
    );                                                                         \
                                                                               \
    /* Composition at an interface with non-pure mixtures */                   \
    makeSpecieNernstModel(Model, Thermo1, Thermo2)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace nernstModels;

    makeNernstModel
    (
        standard,
        gasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeSpecieNernstModel
    (
        standard,
        gasEThermoPhysics,
        gasEThermoPhysics
    );
    makeSpecieNernstModel
    (
        standard,
        constGasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeSpecieNernstModel
    (
        standard,
        constGasEThermoPhysics,
        constGasEThermoPhysics
    );
    makeNernstModel
    (
        standard,
        gasEThermoPhysics,
        icoPoly8EThermoPhysics
    );
    makeSpecieNernstModel
    (
        standard,
        constGasEThermoPhysics,
        icoPoly8EThermoPhysics
    );

    makeNernstModel
    (
        fixedValue,
        gasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeSpecieNernstModel
    (
        fixedValue,
        gasEThermoPhysics,
        gasEThermoPhysics
    );
    makeSpecieNernstModel
    (
        fixedValue,
        constGasEThermoPhysics,
        constFluidEThermoPhysics
    );
    makeSpecieNernstModel
    (
        fixedValue,
        constGasEThermoPhysics,
        constGasEThermoPhysics
    );
    makeNernstModel
    (
        fixedValue,
        gasEThermoPhysics,
        icoPoly8EThermoPhysics
    );
    makeSpecieNernstModel
    (
        fixedValue,
        constGasEThermoPhysics,
        icoPoly8EThermoPhysics
    );
}

// ************************************************************************* //
