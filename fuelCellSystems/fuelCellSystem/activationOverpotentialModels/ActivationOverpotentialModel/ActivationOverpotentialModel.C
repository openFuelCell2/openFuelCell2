/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "ActivationOverpotentialModel.H"

#include "phaseModel.H"
#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "rhoThermo.H"

#include "constants.H"

const Foam::scalar Rgas = Foam::constant::physicoChemical::R.value();
const Foam::scalar F = Foam::constant::physicoChemical::F.value();

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
template<class ThermoType>
const typename Foam::multiComponentMixture<ThermoType>::thermoType&
Foam::ActivationOverpotentialModel<Thermo>::getLocalThermo
(
    const word& speciesName,
    const multiComponentMixture<ThermoType>& globalThermo
) const
{
    return
        globalThermo.getLocalThermo
        (
            globalThermo.species()
            [
                speciesName
            ]
        );
}


template<class Thermo>
template<class ThermoType>
const typename Foam::pureMixture<ThermoType>::thermoType&
Foam::ActivationOverpotentialModel<Thermo>::getLocalThermo
(
    const word& speciesName,
    const pureMixture<ThermoType>& globalThermo
) const
{
    return globalThermo.cellMixture(0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::ActivationOverpotentialModel<Thermo>::ActivationOverpotentialModel
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    activationOverpotentialModel(phase, dict),

    thermo_
    (
        phase.mesh().lookupObject<Thermo>
        (
            IOobject::groupName(basicThermo::dictName, phase.name())
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::ActivationOverpotentialModel<Thermo>::~ActivationOverpotentialModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
