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

#include "ChangJaffe.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ChangJaffe<Thermo>::ChangJaffe
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    ActivationOverpotentialModel<Thermo>(phase, dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ChangJaffe<Thermo>::~ChangJaffe()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::activationOverpotentialModels::ChangJaffe<Thermo>::correct()
{
    //- Create Nernst model
    this->createNernst();

    //- Correct
    this->nernst_->correct();

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    const regionType& fluidPhase = this->region
    (
        word(this->regions_.subDict("fluid").lookup("name"))
    );
    const regionType& electronPhase = this->region
    (
        word(this->regions_.subDict("electron").lookup("name"))
    );
    const regionType& ionPhase = this->region
    (
        word(this->regions_.subDict("ion").lookup("name"))
    );

    //- Source/sink terms for electron and ion fields
    scalarField& SE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>(this->jName)
    );
    scalarField& SI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>(this->jName)
    );
    scalarField& jPrimeE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("jPrime")
    );
    scalarField& jPrimeI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("jPrime")
    );
    scalarField& jSourceE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("jSource")
    );
    scalarField& jSourceI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("jSource")
    );

    //- Potential fields
    const scalarField& phiE = electronPhase.template
        lookupObject<volScalarField>
        (
            this->phiName
        );
    const scalarField& phiI = ionPhase.template
        lookupObject<volScalarField>
        (
            this->phiName
        );

    //- Reference
    scalarField& eta = this->eta_;
    //- Nernst field
    scalarField& nernst = this->nernst_()();
    //- Current density
    scalarField& j = this->j_;

    //- Access reactant
    const scalarField& s = this->phase_;

    //- Mixture mole fraction
    const scalarField W(this->thermo_.W()/1000.0);
    //- Store values for all species

    scalarField coeff(this->phase_.mesh().nCells(), 1.0);
    //- Loop for all relative species
    forAllConstIter(dictionary, this->species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict0 = iter().dict();

            scalar ksi = readScalar(dict0.lookup("ksi"));

            const scalarField& X = this->phase_.X(name);

            //- Change from molar concentration to molar fraction
//            coeff *= Foam::pow(X*rho/W/cRef, ksi);
            coeff *= Foam::pow(X, ksi);
        }
    }

    //- Only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(this->zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    //- electron transfer
    //- cathode side, < 0
    //- anode side, >0
    scalar n = this->nernst_->rxnList()["e"];
    scalar sign = n/mag(n);
    scalar relax = this->relax();
    const scalar j0 = this->j0().value();

    //- The total current: volume averaged
    scalar Rj(0.0);

    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label electronId = electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        //- activation overpotential
        eta[fluidId] = 
        (
            eta[fluidId]*(scalar(1) - relax)
          + (phiE[electronId] - phiI[ionId] - nernst[fluidId])
          * relax
        );

        const scalar jPrime =
            j0*coeff[fluidId]*Foam::pow(s[fluidId], this->gamma_);
        j[fluidId] = jPrime*sign*eta[fluidId];

        //- anode side: SE = Rj, SI = -Rj
        //- cathode side: SE = -Rj, SI = Rj
        SE[electronId] = -sign*j[fluidId];
        SI[ionId] = -SE[electronId];
        jPrimeE[electronId] = jPrime;
        jPrimeI[ionId] = jPrime;
        jSourceE[electronId] =
            SE[electronId]
          + jPrime*(phiI[ionId] + nernst[fluidId] + eta[fluidId]);
        jSourceI[ionId] =
            SI[ionId]
          + jPrime*(phiE[electronId] - nernst[fluidId] - eta[fluidId]);

        Rj += fluidPhase.V()[fluidId] * j[fluidId];
    }

    reduce(Rj, sumOp<scalar>());

    Info << "Total current (A) at " << this->zoneName_
         << ": " << Rj << endl;
}

// ************************************************************************* //
