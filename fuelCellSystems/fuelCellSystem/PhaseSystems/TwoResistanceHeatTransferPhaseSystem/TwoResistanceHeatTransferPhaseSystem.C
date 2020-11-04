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

#include "TwoResistanceHeatTransferPhaseSystem.H"

#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
TwoResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "heatTransfer",
        heatTransferModels_,
        false
    );

    // Check that models have been specified on both sides of the interfaces
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        if (!heatTransferModels_[pair].first().valid())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase1().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
        if (!heatTransferModels_[pair].second().valid())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase2().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
    }

    // Calculate initial Tf-s as if there is no mass transfer
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& T1(phase1.thermo().T());
        const volScalarField& T2(phase2.thermo().T());

        volScalarField H1(heatTransferModels_[pair].first()->K());
        volScalarField H2(heatTransferModels_[pair].second()->K());
        dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

        Tf_.set
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Tf", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (H1*T1 + H2*T2)/max(H1 + H2, HSmall)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~TwoResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.set
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    // Heat transfer with the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const volScalarField& Tf(*Tf_[pair]);

        const Pair<tmp<volScalarField>> Ks
        (
            heatTransferModelIter().first()->K(),
            heatTransferModelIter().second()->K()
        );

        const volScalarField KEff
        (
            Ks.first()()*Ks.second()()
           /max
            (
                Ks.first()() + Ks.second()(),
                dimensionedScalar("small", heatTransferModel::dimK, SMALL)
            )
        );

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();

            const volScalarField& he(phase.thermo().he());
            const volScalarField Cpv(phase.thermo().Cpv());

            *eqns[phase.name()] +=
                Ks[iter.index()]*(Tf - phase.thermo().T())
              + KEff/Cpv*he - fvm::Sp(KEff/Cpv, he);
        }
    }

    // Source term due to mass transfer
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[phase1.name()] += - fvm::Sp(dmdt21, he1) + dmdt21*(K2 - K1);

        *eqns[phase2.name()] -= - fvm::Sp(dmdt12, he2) + dmdt12*(K1 - K2);

        if (this->heatTransferModels_.found(phasePairIter.key()))
        {
            const volScalarField& Tf(*Tf_[pair]);

            *eqns[phase1.name()] +=
                dmdt21*phase1.thermo().he(phase1.thermo().p(), Tf);

            *eqns[phase2.name()] -=
                dmdt12*phase2.thermo().he(phase2.thermo().p(), Tf);
        }
        else
        {
            *eqns[phase1.name()] += dmdt21*he2;

            *eqns[phase2.name()] -= dmdt12*he1;
        }
    }

    return eqnsPtr;
}



template<class BasePhaseSystem>
Foam::autoPtr<Foam::fvScalarMatrix>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::heatTransfer
(
    volScalarField& T,
    const labelList& cellMap
) const
{
    const word& continuous = phaseSystem::continuous();

    autoPtr<fvScalarMatrix> eqnPtr
    (
        new fvScalarMatrix(T, dimEnergy/dimTime)
    );

    const fvMesh& mesh = T.mesh();

    // Heat transfer with the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        if 
        (
            phase1.name() != continuous ||
            phase2.name() != continuous
        )
        {
            continue;
        }

        const volScalarField& Tf(*Tf_[pair]);

        const Pair<tmp<volScalarField>> Ks
        (
            heatTransferModelIter().first()->K(),
            heatTransferModelIter().second()->K()
        );

        const volScalarField KEff
        (
            Ks.first()()*Ks.second()()
           /max
            (
                Ks.first()() + Ks.second()(),
                dimensionedScalar("small", heatTransferModel::dimK, SMALL)
            )
        );

        volScalarField otherT("otherT", 0.0*T);
        otherT.rmap(Tf, cellMap);

        volScalarField KEff0
        (
            IOobject
            (
                "KEff0",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("KEff", heatTransferModel::dimK, 0.0)
        );
        KEff0.rmap(KEff, cellMap);

        if (phase1.name() == continuous)
        {
            volScalarField Ks0("Ks0", KEff0*0.0);
            Ks0.rmap(Ks.first().ref(), cellMap);

            eqnPtr() +=
                Ks0*otherT - fvm::Sp(Ks0, T)
              + KEff0*T - fvm::Sp(KEff0, T);
        }
        else
        {
            volScalarField Ks0("Ks0", KEff0*0.0);
            Ks0.rmap(Ks.second().ref(), cellMap);

            eqnPtr() +=
                Ks0*otherT - fvm::Sp(Ks0, T)
              + KEff0*T - fvm::Sp(KEff0, T);
        }
    }

    // Source term due to mass transfer
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        if 
        (
            phase1.name() != continuous ||
            phase2.name() != continuous
        )
        {
            continue;
        }

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        volScalarField dmdthe0
        (
            IOobject
            (
                "dmdthe",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("dmdthe", dimEnergy/dimTime, 0.0)
        );

        volScalarField dmdtCp0
        (
            IOobject
            (
                "dmdtCp",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("dmdtCp", dimEnergy/dimTime/dimTemperature, 0.0)
        );

        volScalarField dmdtKK0
        (
            IOobject
            (
                "dmdtKK",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("dmdtKK", dimEnergy/dimTemperature, 0.0)
        );

        if (this->heatTransferModels_.found(phasePairIter.key()))
        {
            const volScalarField& Tf(*Tf_[pair]);

            if (phase1.name() == continuous)
            {
                volScalarField dmdthe =
                    dmdt21*phase1.thermo().he(phase1.thermo().p(), Tf);

                dmdthe0.rmap(dmdthe, cellMap);
            }
            else
            {
                volScalarField dmdthe =
                    dmdt12*phase2.thermo().he(phase2.thermo().p(), Tf);

                dmdthe0.rmap(dmdthe, cellMap);
            }
        }
        else
        {
            if (phase1.name() == continuous)
            {
                volScalarField dmdthe = dmdt21*he2;

                dmdthe0.rmap(dmdthe, cellMap);
            }
            else
            {
                volScalarField dmdthe = dmdt12*he1;

                dmdthe0.rmap(dmdthe, cellMap);
            }
        }


        if (phase1.name() == continuous)
        {
            volScalarField dmdtCp = dmdt21*phase1.thermo().Cp().ref();
            volScalarField dmdtKK = dmdt21*(K2 - K1);

            dmdtCp0.rmap(dmdtCp, cellMap);
            dmdtKK0.rmap(dmdtKK, cellMap);

            eqnPtr() += dmdthe0 - fvm::Sp(dmdtCp0, T) + dmdtKK0;
        }
        else
        {
            volScalarField dmdtCp = dmdt12*phase2.thermo().Cp().ref();
            volScalarField dmdtKK = dmdt12*(K1 - K2);

            dmdtCp0.rmap(dmdtCp, cellMap);
            dmdtKK0.rmap(dmdtKK, cellMap);

            eqnPtr() -= dmdthe0 - fvm::Sp(dmdtCp0, T) + dmdtKK0;
        }
    }

    return eqnPtr;
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctThermo()
{
    phaseSystem::correctThermo();

    correctInterfaceThermo();
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& p(phase1.thermo().p());

        const volScalarField& T1(phase1.thermo().T());
        const volScalarField& T2(phase2.thermo().T());

        volScalarField& Tf(*this->Tf_[pair]);

        const volScalarField L
        (
            phase1.thermo().he(p, Tf) - phase2.thermo().he(p, Tf)
        );

        const volScalarField dmdt(this->dmdt(pair));

        volScalarField H1
        (
            this->heatTransferModels_[pair].first()->K()
        );

        volScalarField H2
        (
            this->heatTransferModels_[pair].second()->K()
        );

        // Limit the H[12] to avoid /0
        H1.max(SMALL);
        H2.max(SMALL);

        Tf = (H1*T1 + H2*T2 + dmdt*L)/(H1 + H2);

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.primitiveField())
            << ", mean = " << average(Tf.primitiveField())
            << ", max = " << max(Tf.primitiveField())
            << endl;
    }
}


template<class BasePhaseSystem>
bool Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
