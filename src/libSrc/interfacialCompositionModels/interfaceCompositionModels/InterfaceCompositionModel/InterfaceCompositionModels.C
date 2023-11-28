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

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceCompositionModel.H"
#include "InterfaceCompositionModel.H"
#include "Henry.H"
#include "NonRandomTwoLiquid.H"
#include "Raoult.H"
#include "Saturated.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"

#include "rhoConst.H"
#include "perfectFluid.H"

#include "pureMixture.H"
#include "coefficientMultiComponentMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "forThermo.H"
#include "makeThermo.H"
#include "makeReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSpecieInterfaceCompositionModel(Model, Thermo1, Thermo2)           \
                                                                               \
    /* Composition at an interface with a multi-component mixture */           \
    makeSpecieInterfaceCompositionType                                         \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        coefficientMultiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoReactionThermo,                                        \
        coefficientMultiComponentMixture, Thermo2                                         \
    );

#define makeInterfaceCompositionModel(Model, Thermo1, Thermo2)                 \
                                                                               \
    /* Composition at an interface with a pure mixture */                      \
    makeInterfaceCompositionType                                               \
    (                                                                          \
        Model,                                                                 \
        heRhoThermo, rhoReactionThermo,                                        \
        coefficientMultiComponentMixture, Thermo1,                                        \
        heRhoThermo, rhoThermo,                                                \
        pureMixture, Thermo2                                                   \
    );                                                                         \
                                                                               \
    /* Composition at an interface with non-pure mixtures */                   \
    makeSpecieInterfaceCompositionModel(Model, Thermo1, Thermo2)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace interfaceCompositionModels;

    // Gas-side models
    makeInterfaceCompositionModel
    (
        Saturated,
        gasEThermoPhysics,
        rPolyHThermoPhysics
    );
    makeInterfaceCompositionModel
    (
        Saturated,
        constGasEThermoPhysics,
        rPolyHThermoPhysics
    );

    makeSpecieInterfaceCompositionModel
    (
        NonRandomTwoLiquid,
        gasEThermoPhysics,
        rPolyHThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        NonRandomTwoLiquid,
        constGasEThermoPhysics,
        rPolyHThermoPhysics
    );

    // Liquid-side models
    makeSpecieInterfaceCompositionModel
    (
        Henry,
        rPolyHThermoPhysics,
        gasEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        Henry,
        rPolyHThermoPhysics,
        constGasEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        Raoult,
        rPolyHThermoPhysics,
        gasEThermoPhysics
    );
    makeSpecieInterfaceCompositionModel
    (
        Raoult,
        rPolyHThermoPhysics,
        constGasEThermoPhysics
    );
}

// ************************************************************************* //
