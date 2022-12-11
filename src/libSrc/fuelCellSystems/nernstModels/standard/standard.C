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

#include "standard.H"

#include "constants.H"

const Foam::scalar Rgas = Foam::constant::physicoChemical::R.value();
const Foam::scalar F = Foam::constant::physicoChemical::F.value();

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::nernstModels::standard<Thermo, OtherThermo>::standard
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dictionary& dict
)
:
    NernstModel<Thermo, OtherThermo>(phase1, phase2, dict),
    phase1_(phase1),
    phase2_(phase2)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::nernstModels::standard<Thermo, OtherThermo>::~standard()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::nernstModels::standard<Thermo, OtherThermo>::correct()
{
    this->deltaH_ *= 0.0;
    this->deltaS_ *= 0.0;

    const word water = phaseModel::water;

    scalarField& nernst(*this);
    scalarField Qrxn(this->size(), 1.0);
    scalarField deltaHS(this->size(), 0.0);

    const scalarField& p = this->thermo_.p();
    const scalarField& T = this->thermo_.T();

    scalarField pRef = p/this->pRef_.value();

    forAllConstIter(HashTable<scalar>, this->rxnList_, iter)
    {
        const word& nameI = iter.key();

        if (nameI != "e" && nameI != water)
        {
            label speciesI = this->thermo_.composition().species()[nameI];
            const scalar Wi = this->thermo_.composition().W(speciesI)/1000.0;
            scalar stoiCoeffI = this->rxnList_[nameI];

            forAll(this->deltaH_, cellI)
            {
                this->deltaH_[cellI] +=
                    stoiCoeffI*this->thermo_.composition().Ha(speciesI, p[cellI], T[cellI])*Wi;

                this->deltaS_[cellI] +=
                    stoiCoeffI*this->thermo_.composition().S(speciesI, p[cellI], T[cellI])*Wi;

                deltaHS[cellI] +=
                    stoiCoeffI*this->thermo_.composition().S(speciesI, p[cellI], T[cellI])*Wi*T[cellI];
            }

            const scalarField& X = phase1_.X(nameI);

            Qrxn *= Foam::pow(Foam::max(X, 1.0e-6)*pRef, stoiCoeffI);
        }
    }

    if (this->rxnList_.found(water))
    {
        const scalarField& p = this->otherThermo_.p();
        const scalarField& T = this->otherThermo_.T();

        const typename OtherThermo::thermoType& otherLocalThermo =
            this->getLocalThermo
            (
                water,
                this->otherThermo_
            );

        label speciesI = this->thermo_.composition().species()[water];
        const scalar Wi = this->thermo_.composition().W(speciesI)/1000.0;
        scalar stoiCoeffI = this->rxnList_[water];

        forAll(this->deltaH_, cellI)
        {
            this->deltaH_[cellI] += stoiCoeffI*otherLocalThermo.Ha(p[cellI], T[cellI])*Wi;
            this->deltaS_[cellI] += stoiCoeffI*otherLocalThermo.S(p[cellI], T[cellI])*Wi;
            deltaHS[cellI] += stoiCoeffI*otherLocalThermo.S(p[cellI], T[cellI])*Wi*T[cellI];
        }

        if (phase1_.name() == phase2_.name())
        {
            const scalarField& X = phase1_.X(water);

            Qrxn *= Foam::pow(Foam::max(X, 1.0e-6)*pRef, stoiCoeffI);
        }
    }

    nernst = -(-(this->deltaH_ - deltaHS) - Rgas*T*Foam::log(Qrxn))/this->rxnList_["e"]/F;

    Info<< "Nernst " << this->mesh().name()
        << ": min = " << Foam::min(this->primitiveField())
        << ", mean = " << Foam::average(this->primitiveField())
        << ", max = " << Foam::max(this->primitiveField())
        << endl;
}

// ************************************************************************* //
