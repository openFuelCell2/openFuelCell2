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

#include "fixedValue.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::nernstModels::fixedValue<Thermo, OtherThermo>::fixedValue
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dictionary& dict
)
:
    NernstModel<Thermo, OtherThermo>(phase1, phase2, dict),
    phase1_(phase1),
    phase2_(phase2),
    E0_("E0", this->operator()().dimensions(), dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::nernstModels::fixedValue<Thermo, OtherThermo>::~fixedValue()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::nernstModels::fixedValue<Thermo, OtherThermo>::correct()
{
    this->deltaH() *= 0.0;
    this->deltaS() *= 0.0;

    const word water = phaseModel::water;

    scalarField& nernst = this->operator()();
    scalarField Qrxn(nernst.size(), 1.0);

    const scalarField& p = this->thermo_.p();
    const scalarField& T = this->thermo_.T();

    scalarField pRef = p/this->pRef().value();

    forAllConstIter(HashTable<scalar>, this->rxnList(), iter)
    {
        const word& nameI = iter.key();

        if (nameI != "e" && nameI != water)
        {
            label speciesI = this->thermo_.composition().species()[nameI];
            const scalar Wi = this->thermo_.composition().W(speciesI)/1000.0;
            scalar stoiCoeffI = this->rxnList()[nameI];

            forAll(this->deltaH(), cellI)
            {
                this->deltaH()[cellI] +=
                    stoiCoeffI*this->thermo_.composition().Hc(speciesI)*Wi;

                this->deltaS()[cellI] +=
                    stoiCoeffI*this->thermo_.composition().S(speciesI, p[cellI], T[cellI])*Wi;
            }

            const scalarField& X = phase1_.X(nameI);

            Qrxn *= Foam::pow(Foam::max(X, this->residualY())*pRef, stoiCoeffI);
        }
    }

    if (this->rxnList().found(water))
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
        scalar stoiCoeffI = this->rxnList()[water];

        forAll(this->deltaH(), cellI)
        {
            this->deltaH()[cellI] += stoiCoeffI*otherLocalThermo.Hc()*Wi;
            this->deltaS()[cellI] += stoiCoeffI*otherLocalThermo.S(p[cellI], T[cellI])*Wi;
        }

        if (phase1_.name() == phase2_.name())
        {
            const scalarField& X = phase1_.X(water);

            Qrxn *= Foam::pow(Foam::max(X, this->residualY())*pRef, stoiCoeffI);
        }
    }

    nernst = E0_.value() + Rgas*T*Foam::log(Qrxn)/this->rxnList()["e"]/F;

    Info<< "Nernst " << this->operator()().mesh().name()
        << ": min = " << Foam::min(this->operator()().primitiveField())
        << ", mean = " << Foam::average(this->operator()().primitiveField())
        << ", max = " << Foam::max(this->operator()().primitiveField())
        << endl;
}

// ************************************************************************* //

