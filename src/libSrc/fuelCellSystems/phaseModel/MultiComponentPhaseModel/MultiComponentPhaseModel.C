/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | openFuelCell
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

#include "MultiComponentPhaseModel.H"

#include "phaseSystem.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentPhaseModel<BasePhaseModel>::MultiComponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.mesh().solverDict("Yi")
    ),
    inertIndex_(-1)
{
    const word inertSpecie
    (
        this->thermo_->lookupOrDefault("inertSpecie", word::null)
    );

    if (inertSpecie != word::null)
    {
        inertIndex_ = this->thermo_->composition().species()[inertSpecie];
    }

    PtrList<volScalarField>& Y = this->thermo_->composition().Y();

    //- Mixture mole fraction
    const volScalarField W(this->thermo_->W());

    //- Mole fraction, X
    X_.resize(Y.size());
    iDmdt_.resize(Y.size());

    //- Diffusion coefficients
    diffSp_.resize(Y.size());

    forAll(Y, i)
    {
        const dimensionedScalar Wi
        (
            "W",
            dimMass/dimMoles,
            this->thermo_->composition().W(i)
        );

        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("X", Y[i].member()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                W*Y[i]/Wi
            )
        );

        iDmdt_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdt" + phaseName, Y[i].member()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("iDmdt", dimDensity/dimTime, 0.0)
            )
        );

        if (i != inertIndex_ && this->thermo_->composition().active(i))
        {
            const label j = YActive_.size();
            YActive_.resize(j + 1);
            YActive_.set(j, &Y[i]);

            XActive_.resize(j + 1);
            XActive_.set(j, &X_[i]);
        }

        diffSp_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("diffSp", Y[i].member()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("diff", this->muEff()->dimensions()/dimDensity, SMALL)
            )
        );
    }

    // Diffusivity
    diffModels_.set
    (
        new diffusivityModelList
        (
            this->mesh(),
            this->name()
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentPhaseModel<BasePhaseModel>::~MultiComponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MultiComponentPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    volScalarField Yt
    (
        IOobject
        (
            IOobject::groupName("Yt", this->name()),
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0)
    );

    PtrList<volScalarField>& Yi = YRef();

    forAll(Yi, i)
    {
        if (i != inertIndex_)
        {
            Yt += Yi[i];
        }
    }

    if (inertIndex_ != -1)
    {
        Yi[inertIndex_] = scalar(1) - Yt;
        Yi[inertIndex_].max(0);

        Yt += Yi[inertIndex_];
    }

    //- Normalize
    forAll(Yi, i)
    {
        Yi[i] /= Yt;
    }

    //- Mixture mole fraction
    const volScalarField W(this->thermo_->W());

    //- Update mole fraction X_
    forAll(Yi, i)
    {
        const dimensionedScalar Wi
        (
            "W",
            dimMass/dimMoles,
            this->thermo_->composition().W(i)
        );

        X_[i] = W*Yi[i]/Wi;
    }

    diffModels_->correct(X_, diffSp_, inertIndex_);
}


template<class BasePhaseModel>
bool Foam::MultiComponentPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    const volScalarField& alpha = *this;
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());
    const volScalarField& rho = this->thermo().rho();

    label speciesIndex = this->thermo_->composition().species()[Yi.member()];

    return
    (
        fvm::ddt(alpha, rho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")

      - fvm::laplacian
        (
            fvc::interpolate(alpha*rho)
           *fvc::interpolate(diffSp_[speciesIndex]),
            Yi,
            "laplacian(diff,Yi)"
        )
      ==
        alpha*this->R(Yi)

      + fvc::ddt(residualAlpha_*rho, Yi)
      - fvm::ddt(residualAlpha_*rho, Yi)
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::iDmdt() const
{
    return iDmdt_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::X() const
{
    return X_;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    return this->thermo_->composition().Y(name);
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentPhaseModel<BasePhaseModel>::iDmdt(const word& name) const
{
    return iDmdt_[this->thermo_->composition().species()[name]];
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentPhaseModel<BasePhaseModel>::X(const word& name) const
{
    return X_[this->thermo_->composition().species()[name]];
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YRef()
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::iDmdtRef()
{
    return iDmdt_;
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::XRef()
{
    return X_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YActive() const
{
    return YActive_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::XActive() const
{
    return XActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YActiveRef()
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::XActiveRef()
{
    return XActive_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel>::dmdt() const
{
    tmp<volScalarField> dmdt
    (
        volScalarField::New
        (
            IOobject::groupName("dmdt0", this->name()),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );

    scalarField& dmdt0 = dmdt.ref();

    for (auto& iDmdt : iDmdt_)
    {
        dmdt0 += iDmdt;
    }

    return dmdt;
}

// ************************************************************************* //
