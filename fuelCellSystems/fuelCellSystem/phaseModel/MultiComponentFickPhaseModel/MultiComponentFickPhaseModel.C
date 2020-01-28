/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "MultiComponentFickPhaseModel.H"

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
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::MultiComponentFickPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    diff_(fluid.mesh().nCells(), 0.0),
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
                    IOobject::NO_READ,
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

    //- Pointerlist for all diffusivity models
    diffModels_.setSize(1 + fluid.porousZone().size());

    // Diffusivity: whole zone
    // ---------------------------
    label znId = this->mesh().cellZones().findZoneID(this->mesh().name());
    const labelList& cells = this->mesh().cellZones()[znId];
    const dictionary& diffDict = fluid.subDict(this->name()).subDict("diffusivity");

    diffModels_.set
    (
        0,
        new autoPtr<diffusivityModels::diffusivityModel>
        (
            diffusivityModels::diffusivityModel::New
            (
                this->mesh(),
                diff_,
                cells, 
                diffDict
            )
        )
    );

    // porous zones
    forAll(fluid.porousZone(), iz)
    {
        label znID =
              this->mesh().cellZones().findZoneID(fluid.porousZone()[iz].zoneName());

        const labelList& cells = this->mesh().cellZones()[znID];
        const dictionary& porousDiffDict = fluid.porousZone()[iz].dict().subDict("diffusivity");

        diffModels_.set
        (
            iz + 1,
            new autoPtr<diffusivityModels::diffusivityModel>
            (
                diffusivityModels::diffusivityModel::New
                (
                    this->mesh(),
                    diff_,
                    cells, 
                    porousDiffDict
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::~MultiComponentFickPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MultiComponentFickPhaseModel<BasePhaseModel>::correctThermo()
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

    //- Diffusivity model update
    forAll(Yi, i)
    {
        //diffSp[i] is diffusivity field for specie i, to be used in YEqn
        scalarField& diffSpIn = diffSp_[i];
        diffSpIn = 0;

        if(i != inertIndex_)
        {
            Info<< nl << "species " << Yi[i].member() << nl;

            forAll(diffModels_, j)
            {
                if(diffModels_[j]->isFixed())
                {
                    diffModels_[j]->evaluate();

                    diffModels_[j]->writeData();

                    // copy calculated zone diffusivity from diffModel to diffSp[i]
                    forAll(diffModels_[j]->cells(), k)
                    {
                        label iCell = diffModels_[j]->cells()[k];
                        diffSpIn[iCell] = diffModels_[j]->diff()[iCell];
                    }
                }
                else if(!diffModels_[j]->isBinary())
                {
                    FatalErrorIn("ERROR: multicomponent diffusivity")
                        << "requires fixed or binary model"
                        << exit(FatalError);
                }
                else
                {
                    // species diffusivities in mixture from binary diffusivities
                    // ----------------------------------------------------------
                    // D_{i} = (1-x[i])/sum_{l!=i}(x[l]/D[i,l])

                    // pairwise binary diff calculation and accumulation
                    // -------------------------------------------
                    //initialize sum(x[l]/D[il]
                    scalarField sum(diffSpIn.size(), 0);

                    forAll(Yi, l)
                    {
                        if (l != i)
                        {
                            diffModels_[j]->setSpecies
                            (
                                Yi[i].member(),
                                Yi[l].member()
                            );
                            diffModels_[j]->evaluate();
                            diffModels_[j]->writeData();

                            forAll(diffModels_[j]->cells(), k)
                            {
                                label iCell = diffModels_[j]->cells()[k];
                                if(diffModels_[j]->diff()[iCell] != 0) //hkr: surely this must be!
                                {
                                    sum[iCell] += X_[l][iCell]/diffModels_[j]->diff()[iCell];
                                }
                            }
                        }
                    }

                    // diffSpIn[zone] <-- (1-x[a])/sum
                    forAll(diffModels_[j]->cells(), k)
                    {
                        label iCell = diffModels_[j]->cells()[k];
                        if (sum[iCell] != 0)
                        {
                            diffSpIn[iCell] = (1 - X_[i][iCell])/sum[iCell];
                        }
                    }
                } //isBinary
            } //m
            diffSp_[i].correctBoundaryConditions();
        } //!inert
    }
}


template<class BasePhaseModel>
bool Foam::MultiComponentFickPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
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
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::iDmdt() const
{
    return iDmdt_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::X() const
{
    return X_;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    return this->thermo_->composition().Y(name);
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::iDmdt(const word& name) const
{
    return iDmdt_[this->thermo_->composition().species()[name]];
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::X(const word& name) const
{
    return X_[this->thermo_->composition().species()[name]];
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::YRef()
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::iDmdtRef()
{
    return iDmdt_;
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::XRef()
{
    return X_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::YActive() const
{
    return YActive_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::XActive() const
{
    return XActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::YActiveRef()
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::XActiveRef()
{
    return XActive_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MultiComponentFickPhaseModel<BasePhaseModel>::dmdt() const
{
    tmp<volScalarField> dmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dmdt0", this->name()),
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    scalarField& dmdt0 = dmdt.ref();

    forAll(Y(), i)
    {
        volScalarField& Yi = const_cast<volScalarField&>(Y()[i]);

        dmdt0 += iDmdt_[i];
    }

    return dmdt;
}

// ************************************************************************* //
