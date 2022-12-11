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

#include "NernstModel.H"

#include "Pair.H"
#include "phaseSystem.H"
#include "phasePair.H"
#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "rhoThermo.H"

#include "constants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::multiComponentMixture<ThermoType>::thermoType&
Foam::NernstModel<Thermo, OtherThermo>::getLocalThermo
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


template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::pureMixture<ThermoType>::thermoType&
Foam::NernstModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const pureMixture<ThermoType>& globalThermo
) const
{
    return globalThermo.cellMixture(0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::NernstModel<Thermo, OtherThermo>::NernstModel
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dictionary& dict
)
:
    nernstModel(phase1, phase2, dict),

    thermo_
    (
        phase1.mesh().lookupObject<Thermo>
        (
            IOobject::groupName(basicThermo::dictName, phase1.name())
        )
    ),
    otherThermo_
    (
        phase2.mesh().lookupObject<OtherThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase2.name())
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::NernstModel<Thermo, OtherThermo>::~NernstModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
