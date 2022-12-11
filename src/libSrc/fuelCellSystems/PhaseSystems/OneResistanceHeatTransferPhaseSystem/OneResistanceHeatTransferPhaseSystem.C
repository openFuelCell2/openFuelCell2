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

#include "OneResistanceHeatTransferPhaseSystem.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
OneResistanceHeatTransferPhaseSystem
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~OneResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
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

    // Heat transfer across the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const volScalarField K(heatTransferModelIter()->K());

        const phasePair& pair(this->phasePairs_[heatTransferModelIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField& he(phase.thermo().he());
            volScalarField Cpv(phase.thermo().Cpv());

            *eqns[phase.name()] +=
                K*(otherPhase.thermo().T() - phase.thermo().T() + he/Cpv)
              - fvm::Sp(K/Cpv, he);
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

        // Note that the phase heEqn contains a continuity error term, which
        // implicitly adds a mass transfer term of fvm::Sp(dmdt, he). These
        // additions do not include this term.

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[phase1.name()] +=
            dmdt21*he2 - fvm::Sp(dmdt21, he1) + dmdt21*(K2 - K1);

        *eqns[phase2.name()] -=
            dmdt12*he1 - fvm::Sp(dmdt12, he2) + dmdt12*(K1 - K2);
    }

    return eqnsPtr;
}



template<class BasePhaseSystem>
Foam::autoPtr<Foam::fvScalarMatrix>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer
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

    // Heat transfer across the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const volScalarField K(heatTransferModelIter()->K());

        const phasePair& pair(this->phasePairs_[heatTransferModelIter.key()]);

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        if 
        (
            phase1.name() != continuous &&
            phase2.name() != continuous
        )
        {
            continue;
        }

        volScalarField KT
        (
            IOobject
            (
                "KT",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("k", K.dimensions(), 0.0)
        );

        volScalarField otherT("otherT", 0.0*T);

        if (phase1.name() == continuous)
        {
            scalarField T0 = phase2.thermo().T();

            otherT.rmap(T0, cellMap);
        }
        else
        {
            scalarField T0 = phase1.thermo().T();

            otherT.rmap(T0, cellMap);
        }

        KT.rmap(K, cellMap);

        eqnPtr() += KT*otherT - fvm::Sp(KT, T);
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
            phase1.name() != continuous &&
            phase2.name() != continuous
        )
        {
            continue;
        }

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Note that the phase heEqn contains a continuity error term, which
        // implicitly adds a mass transfer term of fvm::Sp(dmdt, he). These
        // additions do not include this term.

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
            dimensionedScalar("dmdthe", dimPower/dimVol, 0.0)
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
            dimensionedScalar("dmdtCp", dimPower/dimTemperature, 0.0)
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
            dimensionedScalar("dmdtKK", dimPower/dimVol, 0.0)
        );

        if (phase1.name() == continuous)
        {
            volScalarField dmdtCp = dmdt21*phase1.thermo().Cp().ref();
            volScalarField dmdtKK = dmdt21*(K2 - K1);
            volScalarField dmdthe = dmdt21*he2;

            dmdthe0.rmap(dmdthe, cellMap);
            dmdtCp0.rmap(dmdtCp, cellMap);
            dmdtKK0.rmap(dmdtKK, cellMap);

            eqnPtr() += dmdthe0 - fvm::Sp(dmdtCp0, T) + dmdtKK0;
        }
        else
        {
            volScalarField dmdtCp = dmdt12*phase2.thermo().Cp().ref();
            volScalarField dmdtKK = dmdt12*(K1 - K2);
            volScalarField dmdthe = dmdt12*he1;

            dmdthe0.rmap(dmdthe, cellMap);
            dmdtCp0.rmap(dmdtCp, cellMap);
            dmdtKK0.rmap(dmdtKK, cellMap);

            eqnPtr() -= dmdthe0 - fvm::Sp(dmdtCp0, T) + dmdtKK0;
        }
    }

    return eqnPtr;
}


template<class BasePhaseSystem>
bool Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
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
