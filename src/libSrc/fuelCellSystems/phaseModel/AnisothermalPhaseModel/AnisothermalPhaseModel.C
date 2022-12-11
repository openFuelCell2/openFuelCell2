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

#include "AnisothermalPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::filterPressureWork
(
    const tmp<volScalarField>& pressureWork
) const
{
    const volScalarField& alpha = *this;

    scalar pressureWorkAlphaLimit =
        this->thermo_->lookupOrDefault("pressureWorkAlphaLimit", 0.0);

    if (pressureWorkAlphaLimit > 0)
    {
        return
        (
            max(alpha - pressureWorkAlphaLimit, scalar(0))
           /max(alpha - pressureWorkAlphaLimit, pressureWorkAlphaLimit)
        )*pressureWork;
    }
    else
    {
        return pressureWork;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::AnisothermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::AnisothermalPhaseModel<BasePhaseModel>::~AnisothermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::AnisothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    this->thermo_->correct();
}


template<class BasePhaseModel>
bool Foam::AnisothermalPhaseModel<BasePhaseModel>::isothermal() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;

    const volVectorField U(this->U());
    const surfaceScalarField alphaPhi(this->alphaPhi());
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());

    const volScalarField contErr(this->continuityError());
    const volScalarField K(this->K());

    volScalarField& he = this->thermo_->he();

    //- Porosity
    //- Non-dimension, created via alpha
    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar("por", dimless, 1.0)
    );

    if (!this->fluid().porousZone().empty())
    {
        forAll(this->fluid().porousZone(), iz)
        {
            label znId = this->mesh().
                cellZones().findZoneID(this->fluid().porousZone()[iz].zoneName());

            scalar por = this->fluid().porousZone()[iz].porosity();

            labelList znCells(this->mesh().cellZones()[znId]);

            forAll(znCells, cellI)
            {
                label cellId = znCells[cellI];

                porosity[cellId] *= por;
            }
        }
    }

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha, this->rho(), he)
      + fvm::div(alphaRhoPhi, he)
      + fvm::SuSp(-contErr, he)

      + fvc::ddt(alpha, this->rho(), K) + fvc::div(alphaRhoPhi, K)
      - contErr*K

      - fvm::laplacian
        (
            fvc::interpolate(alpha*porosity)
           *fvc::interpolate(this->alphaEff()),
            he
        )
//     ==
//        alpha*this->Qdot()
    );

    // Add the appropriate pressure-work term
//    if (he.name() == this->thermo_->phasePropertyName("e"))
//    {
//        tEEqn.ref() += filterPressureWork
//        (
//            fvc::div(fvc::absolute(alphaPhi, alpha, U), this->thermo().p())
//          + this->thermo().p()*fvc::ddt(alpha)
//        );
//    }
//    else if (this->thermo_->dpdt())
//    {
//        tEEqn.ref() -= filterPressureWork(alpha*this->fluid().dpdt());
//    }

    return tEEqn;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heQdot()
{
    const volScalarField& alpha = *this;

    const volVectorField U(this->U());
    const surfaceScalarField alphaPhi(this->alphaPhi());
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());

    const volScalarField contErr(this->continuityError());
    const volScalarField K(this->K());

    volScalarField& he = this->thermo_->he();

    tmp<volScalarField> qdot
    (
      - fvc::ddt(alpha, this->rho(), K) - fvc::div(alphaRhoPhi, K)
      + contErr*K
//
//      + alpha*this->Qdot()
    );

    // Add the appropriate pressure-work term
    if (he.name() == this->thermo_->phasePropertyName("e"))
    {
        qdot.ref() -= filterPressureWork
        (
            fvc::div(fvc::absolute(alphaPhi, alpha, U), this->thermo().p())
          + this->thermo().p()*fvc::ddt(alpha)
        );
    }
    else if (this->thermo_->dpdt())
    {
        qdot.ref() += filterPressureWork(alpha*this->fluid().dpdt());
    }

    return qdot;
}

// ************************************************************************* //
