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

#include "electroChemicalReaction.H"

#include "phaseSystem.H"
#include "activationOverpotentialModel.H"
#include "dissolvedModel.H"

#include "constants.H"

const Foam::dimensionedScalar Rgas = Foam::constant::physicoChemical::R;
const Foam::dimensionedScalar dimF = Foam::constant::physicoChemical::F;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class ReactionThermo>
Foam::combustionModels::electroChemicalReaction<ReactionThermo>::electroChemicalReaction
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& electroChemicalReactionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    thermo_(thermo),
    saturation_(saturationModel::New(this->subDict("saturation"), this->mesh())),
    dissolved_(this->template lookupOrDefault<Switch>("dissolved", true))
{
    //- Get current phase Model
    const phaseModel& phase = this->mesh().template
        lookupObject<phaseModel>
        (
            this->thermo_.phasePropertyName("alpha")
        );

    //- Activation overpotential model
    eta_ = activationOverpotentialModel::New(phase, this->subDict("activationOverpotentialModel"));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::electroChemicalReaction<ReactionThermo>::~electroChemicalReaction()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class ReactionThermo>
void Foam::combustionModels::electroChemicalReaction<ReactionThermo>::correct()
{
    //- Activation overpotential model update
    eta_->correct();

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    const regionType& fluidPhase = eta_->region
    (
        word(eta_->regions().subDict("fluid").lookup("name"))
    );
    const regionType& ionPhase = eta_->region
    (
        word(eta_->regions().subDict("ion").lookup("name"))
    );

     // Get the phase System
     const phaseSystem& phaseSys = this->mesh().template
     lookupObject<phaseSystem>(phaseSystem::propertiesName);

    //- Get the present phase Model from fluid phase
    const phaseModel& phase = this->mesh().template
        lookupObject<phaseModel>(thermo_.phasePropertyName("alpha"));

    word water(phaseModel::water);

    //- Lable of water (H2O)
    const label specieI =
        thermo_.composition().species()[water];

    //- In case no water exists
    if (specieI == -1)
    {
        return;
    }

    //- Get the specieStoichCoeff
    dimensionedScalar specieStoichCoeff
    (
        "stoichCoeff",
        dimless,
        eta_->nernst().rxnList().found(water)
      ? eta_->nernst().rxnList()[water]/mag(eta_->nernst().rxnList()["e"])
      : 0.0
    );

    //- kg/kmol -> kg/mol
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        thermo_.composition().W(specieI)/1000
    );

    //- water production
    volScalarField wSpecie
    (
        eta_->j()*Wi*specieStoichCoeff/dimF
    );

    // TODO: A better way to handle the dissolved water transfer.

    //- Dissolved water model
    dissolvedModel& dW = const_cast<dissolvedModel&>
    (
        ionPhase.template
        lookupObject<dissolvedModel>(dissolvedModel::modelName)
    );

    //- Water activity and water source/sink
    scalarField& act = const_cast<volScalarField&>(dW.act());
    scalarField& dmdt = const_cast<volScalarField&>(dW.dmdt());

    //- Water mole fraction
    const scalarField& xH2O = phase.X(water);

    //- Activity
    scalarField act0 =
        xH2O
      * this->thermo_.p()
      / saturation_->pSat(thermo_.T()).ref()
      + 2.*(scalar(1) - phase.primitiveField());

    //- Only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(eta_->zoneName());
    const labelList& cells = fluidPhase.cellZones()[znId];

    //- Update water activity
    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        //- Water activity: act = p(H2O)/pSat + 2*sat
        act[ionId] = act0[fluidId];
    }

    //- Update the source/sink term dmdt
    dW.update(eta_->zoneName());

    //- Update the water production in phase model
    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        if (dissolved_)
        {
            dmdt[ionId] += wSpecie[fluidId]/Wi.value();
        }

        if (phaseSys.isSinglePhase() || !eta_->phaseChange())
        {
            //- Water is transferred between current phase and dissolved phase
            volScalarField& iDmdtWater = const_cast<volScalarField&>
                (phase.iDmdt(water));

            iDmdtWater[fluidId] = wSpecie[fluidId] - dmdt[ionId]*Wi.value();
        }
        else
        {
            //- Get the name of the other phase
            const word name1 = Pair<word>
            (
                phaseSys.phases()[0].name(),
                phaseSys.phases()[1].name()
            ).other(phase.name());

            //- Water is transferred between the other phase and dissolved phase
            scalarField& iDmdtWater = const_cast<volScalarField&>
                (phaseSys.phases()[name1].iDmdt(water));

            iDmdtWater[fluidId] = wSpecie[fluidId] - dmdt[ionId]*Wi.value();
        }
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::electroChemicalReaction<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    //- Get current phase Model
    const phaseModel& phase = this->mesh().template
        lookupObject<phaseModel>
        (
            thermo_.phasePropertyName("alpha")
        );

    word water(phaseModel::water);

    if (Y.member() == water)
    {
        return fvm::Sp(phase.iDmdt(water)/(Y + SMALL), Y);      // TODO: is this a good approach.
    }

    const label specieI =
        thermo_.composition().species()[Y.member()];

    //- kg/kmol -> kg/mol
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        thermo_.composition().W(specieI)/1000
    );

    dimensionedScalar specieStoichCoeff
    (
        "stoichCoeff",
        dimless,
        eta_->nernst().rxnList().found(Y.member())
      ? eta_->nernst().rxnList()[Y.member()]/mag(eta_->nernst().rxnList()["e"])
      : 0.0
    );

    volScalarField wSpecie
    (
        eta_->j()*Wi*specieStoichCoeff/dimF
    );

    volScalarField& iDmdt = const_cast<volScalarField&>(phase.iDmdt(Y.member()));

    iDmdt = wSpecie;

    return fvm::Sp(wSpecie/(Y + SMALL), Y);     // TODO: is this a good approach.
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::electroChemicalReaction<ReactionThermo>::Qdot() const
{
    return eta_->Qdot();
}


template<class ReactionThermo>
bool Foam::combustionModels::electroChemicalReaction<ReactionThermo>::read()
{
    if (ThermoCombustion<ReactionThermo>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
