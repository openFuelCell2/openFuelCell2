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

#include "activationOverpotentialModel.H"
#include "phaseSystem.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Static members  * * * * * * * * * * * * * * //

const Foam::scalar F = Foam::constant::physicoChemical::F.value();
const Foam::word Foam::activationOverpotentialModel::phiName("phi");
const Foam::word Foam::activationOverpotentialModel::jName("J");

namespace Foam
{
    defineTypeNameAndDebug(activationOverpotentialModel, 0);
    defineRunTimeSelectionTable(activationOverpotentialModel, dictionary);
}

// * * * * * * * * * * * * * * * * Protected function  * * * * * * * * * * * * * * //

void Foam::activationOverpotentialModel::createNernst()
{
    //- Create nernst model
    if (!nernst_.valid())
    {
        const phaseSystem& phaseSys = phase_.mesh().template
            lookupObject<phaseSystem>(phaseSystem::propertiesName);

        if (phaseChange_)
        {
            //- Get the name of the other phase
            const word name = Pair<word>
            (
                phaseSys.phases()[0].name(),
                phaseSys.phases()[1].name()
            ).other(phase_.name());

            const phaseModel& phase1 = phaseSys.phases()[name];
            nernst_ = nernstModel::New
            (
                phase_, phase1, dict_.subDict("nernstModel")
            );
        }
        else
        {
            nernst_ = nernstModel::New
            (
                phase_, phase_, dict_.subDict("nernstModel")
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::activationOverpotentialModel::activationOverpotentialModel
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    phase_(phase),
    dict_(dict),
    eta_
    (
        IOobject
        (
            IOobject::groupName("eta", phase.mesh().name()),
            phase.mesh().time().timeName(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar
        (
            "eta",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),
    phaseChange_(dict_.lookupOrDefault<Switch>("phaseChange", true)),
    zoneName_(dict_.lookup("zoneName")),
    regions_(dict_.subDict("regions")),
    species_(dict_.subDict("species")),
    phiNames_(dict_.lookup("phiNames")),
    relax_(dict_.lookupOrDefault<scalar>("relax", 1.0)),
    alpha_(dict_.lookupOrDefault<scalar>("alpha", 0.5)),
    gamma_(dict_.lookupOrDefault<scalar>("gamma", 1.0)),
    j0_("j0", dimCurrent/dimVolume, dict_),

    j_
    (
        IOobject
        (
            "j",
            phase_.mesh().time().timeName(),
            phase_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.mesh(),
        dimensionedScalar
        (
            "j",
            dimCurrent/dimVolume,
            0.0
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationOverpotentialModel::~activationOverpotentialModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::activationOverpotentialModel::Qdot() const
{
    tmp<volScalarField> q
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("q", phase_.name()),
                eta_.mesh().time().timeName(),
                eta_.mesh()
            ),
            eta_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (!nernst_.valid())
    {
        return q;
    }

    //- get the sub models
    const regionType& fluidPhase = region
    (
        word(regions_.subDict("fluid").lookup("name"))
    );

    const scalarField& j = j_;
    const scalarField& T = phase_.thermo().T();

    const scalarField& eta = eta_;

    //- only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    scalar n = nernst_->rxnList()["e"];
    scalar sign = n/mag(n);

    forAll(cells, cellI)
    {
        label fluidId = cells[cellI];

        q.ref()[fluidId] =
            j[fluidId]*(sign*eta[fluidId] - T[fluidId]*nernst_->deltaS()[fluidId]
            /mag(nernst_->rxnList()["e"])/F);
    }

    return q;
}

// ************************************************************************* //
