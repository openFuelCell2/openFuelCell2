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

#include "phaseSystem.H"
#include "driftFluxSystem.H"
#include "OneResistanceHeatTransferPhaseSystem.H"
#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "PhaseTransferPhaseSystem.H"
#include "PopulationBalancePhaseSystem.H"
#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        PhaseTransferPhaseSystem
        <
            OneResistanceHeatTransferPhaseSystem<driftFluxSystem>
        >
        basicDriftFluxSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        basicDriftFluxSystem,
        dictionary,
        basicDriftFluxSystem
    );

    typedef
        InterfaceCompositionPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem<driftFluxSystem>
            >
        >
        interfaceCompositionPhaseChangeDriftFluxSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        interfaceCompositionPhaseChangeDriftFluxSystem,
        dictionary,
        interfaceCompositionPhaseChangeDriftFluxSystem
    );

    typedef
        ThermalPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem<driftFluxSystem>
            >
        >
        thermalPhaseChangeDriftFluxSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        thermalPhaseChangeDriftFluxSystem,
        dictionary,
        thermalPhaseChangeDriftFluxSystem
    );

    typedef
        PopulationBalancePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem<driftFluxSystem>
            >
        >
        populationBalanceDriftFluxSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceDriftFluxSystem,
        dictionary,
        populationBalanceDriftFluxSystem
    );

    typedef
        ThermalPhaseChangePhaseSystem
        <
            PopulationBalancePhaseSystem
            <
                PhaseTransferPhaseSystem
                <
                    TwoResistanceHeatTransferPhaseSystem<driftFluxSystem>
                >
            >
        >
        thermalPhaseChangePopulationBalanceDriftFluxSystem;

        addNamedToRunTimeSelectionTable
        (
            phaseSystem,
            thermalPhaseChangePopulationBalanceDriftFluxSystem,
            dictionary,
            thermalPhaseChangePopulationBalanceDriftFluxSystem
        );
}


// ************************************************************************* //
