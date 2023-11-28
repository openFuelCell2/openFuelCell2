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

#include "Tafel.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::Tafel<Thermo>::Tafel
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    ActivationOverpotentialModel<Thermo>(phase, dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::Tafel<Thermo>::~Tafel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::activationOverpotentialModels::Tafel<Thermo>::correct()
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
    //- See names in regions/electronIon/electronIon.C
    scalarField& SE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("J")
    );
    scalarField& SI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("J")
    );

    //- Potential fields
    const scalarField& phiE = electronPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["electron"])
        );
    const scalarField& phiI = ionPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["ion"])
        );

    //- Reference
    scalarField& eta = this->eta_;
    //- Nernst field
    scalarField& nernst = this->nernst_()();
    //- Current density
    scalarField& j = this->j_;

    //- Access temperature and reactant
    const scalarField& T = this->thermo_.T();
    const scalarField& s = this->phase_;

    //- Mixture mole fraction
    const scalarField W(this->phase_.thermo().W()/1000.0);
    //- Mixture density
    const scalarField& rho = this->phase_.thermo().rho();
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
            scalar cRef = readScalar(dict0.lookup("cRef"));

            const scalarField& X = this->phase_.X(name);

            coeff *= Foam::pow(X*rho/W/cRef, ksi);
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
            eta[fluidId]*(scalar(1) - this->relax_)
          + (phiE[electronId] - phiI[ionId] - nernst[fluidId])
          * this->relax_
        );

        //- Buttler-volmer relation
        j[fluidId] = Foam::max
        (
            this->j0_.value()*
            coeff[fluidId]*
            Foam::pow(s[fluidId], this->gamma_)*
            (
                Foam::exp(n*this->alpha_*F*eta[fluidId]/Rgas/T[fluidId])
            )
            ,
            scalar(0)
        );

        //- anode side: SE = Rj, SI = -Rj
        //- cathode side: SE = -Rj, SI = Rj
        SE[electronId] = -sign*j[fluidId];
        SI[ionId] = -SE[electronId];

        Rj += fluidPhase.V()[fluidId] * j[fluidId];
    }
}

// ************************************************************************* //
