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

#include "fvCFD.H"

#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "virtualMassModel.H"

#include "MULES.H"
#include "subCycle.H"
#include "UniformField.H"

#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseSystem, 0);
}

// * * * * * * * * * * * * * * *Private Functions * * * * * * * * * * * * * //
void Foam::twoPhaseSystem::solveAlpha(phaseModel& phase1)
{
    const Time& runTime = mesh_.time();

    // Define the other phase as phase2
    phaseModel& phase2 = const_cast<phaseModel&>(otherPhase(phase1));

    volScalarField& alpha1 = phase1;
    volScalarField& alpha2 = phase2;

    const dictionary& alphaControls = mesh_.solverDict(alpha1.name());

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    bool LTS = fv::localEulerDdt::enabled(mesh_);

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    const surfaceScalarField& phi1 = phase1.phi();
    const surfaceScalarField& phi2 = phase2.phi();

    // Construct the dilatation rate source term
    tmp<volScalarField::Internal> tdgdt;

    if (phase1.divU().valid() && phase2.divU().valid())
    {
        tdgdt =
        (
            alpha2()
           *phase1.divU()()()
          - alpha1()
           *phase2.divU()()()
        );
    }
    else if (phase1.divU().valid())
    {
        tdgdt =
        (
            alpha2()
           *phase1.divU()()()
        );
    }
    else if (phase2.divU().valid())
    {
        tdgdt =
        (
          - alpha1()
           *phase2.divU()()()
        );
    }

    alpha1.correctBoundaryConditions();

    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alphaDbyA;
    if (DByAfs().found(phase1.name()) && DByAfs().found(phase2.name()))
    {
        surfaceScalarField DbyA
        (
            *DByAfs()[phase1.name()] + *DByAfs()[phase2.name()]
        );

        alphaDbyA =
            fvc::interpolate(max(alpha1, scalar(0)))
           *fvc::interpolate(max(alpha2, scalar(0)))
           *DbyA;

        phir += DbyA*fvc::snGrad(alpha1, "bounded")*mesh_.magSf();
    }
    Info << "Solving for alpha: " << endl;
    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dimless/dimTime, 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        if (tdgdt.valid())
        {
            scalarField& dgdt = tdgdt.ref();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] > 0.0)
                {
                    Sp[celli] -= dgdt[celli]/max(1 - alpha1[celli], 1e-4);
                    Su[celli] += dgdt[celli]/max(1 - alpha1[celli], 1e-4);
                }
                else if (dgdt[celli] < 0.0)
                {
                    Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
                }
            }
        }

        surfaceScalarField alphaPhi1
        (
            fvc::flux
            (
                phi_,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        phase1.correctInflowOutflow(alphaPhi1);

        if (nAlphaSubCycles > 1)
        {
            tmp<volScalarField> trSubDeltaT;

            if (LTS)
            {
                trSubDeltaT =
                    fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
            }

            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhi10(alphaPhi1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi_,
                    alphaPhi10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    UniformField<scalar>(phase1.alphaMax()),
                    zeroField()
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1.alphaPhiRef() = alphaPhi10;
                }
                else
                {
                    phase1.alphaPhiRef() += alphaPhi10;
                }
            }

            phase1.alphaPhiRef() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhi1,
                Sp,
                Su,
                UniformField<scalar>(phase1.alphaMax()),
                zeroField()
            );

            phase1.alphaPhiRef() = alphaPhi1;
        }

        if (alphaDbyA.valid())
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alphaDbyA(), alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1.alphaPhiRef() += alpha1Eqn.flux();
        }

        phase1.alphaRhoPhiRef() =
            fvc::interpolate(phase1.rho())*phase1.alphaPhi();

        phase2.alphaPhiRef() = phi_ - phase1.alphaPhi();
        phase2.correctInflowOutflow(phase2.alphaPhiRef());
        phase2.alphaRhoPhiRef() =
            fvc::interpolate(phase2.rho())*phase2.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(alpha." << alpha1.group() << ") = " << min(alpha1).value()
            << "  Max(alpha." << alpha1.group() << ") = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.clip(phase1.residualAlpha().value(), 1);

        // Update the phase-fraction of the other phase
        alpha2 = scalar(1) - alpha1;
        alpha2.clip(phase2.residualAlpha().value(), 1);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),
    phase1_(phaseModels_[0]),
    phase2_(phaseModels_[1]),
    pc_
    (
        IOobject
        (
            "pc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "pc",
            p_rgh_.dimensions(),
            tensor::zero
        )
    ),

    dpcds_
    (
        IOobject
        (
            "dpcds",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pc_
    )
{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);

    volScalarField& alpha1 = phase1_;
    mesh.setFluxRequired(alpha1.name());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::sigma() const
{
    return sigma
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Kd() const
{
    return Kd
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Kdf() const
{
    return Kdf
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Vm() const
{
    return Vm
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::phiMagMax() const
{
    return max(mag(phase1_.phi()), mag(phase2_.phi()));
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::phirMag() const
{
    return mag(phase1_.phi() - phase2_.phi());
}


void Foam::twoPhaseSystem::solve()
{
    Info << "\nSolving for two phase flow:" << endl;

    fvMesh& mesh = const_cast<fvMesh&>(mesh_);

    const Time& runTime = mesh.time();

    #include "createFields.H"
    #include "createFieldRefs.H"

    Switch faceMomentum
    (
        pimple_.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    word activePhase
    (
        pimple_.dict().lookupOrDefault<word>("activePhase", phase1_.name())
    );

    #include "pUf/createRDeltaTf.H"

    if (LTS && faceMomentum)
    {
        #include "setRDeltaTf.H"
    }

    int nEnergyCorrectors
    (
        pimple_.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
    );

    Switch implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).lookupOrDefault<Switch>
        (
            "implicitPhasePressure", false
        )
    );

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
        solveAlpha(phaseModels_[activePhase]);

        correct();

        #include "YEqns.H"

        if (faceMomentum)
        {
            #include "pUf/UEqns.H"
            #include "EEqn.H"
            #include "pUf/pEqn.H"
        }
        else
        {
            #include "pU/UEqns.H"
            #include "EEqn.H"
            #include "pU/pEqn.H"
        }

        correctKinematics();

        if (pimple_.turbCorr())
        {
            correctTurbulence();
        }
    }
}


bool Foam::twoPhaseSystem::isSinglePhase() const
{
    return false;
}

// ************************************************************************* //
