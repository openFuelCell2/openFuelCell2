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
#include "driftFluxSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "fixedValueFvsPatchFields.H"
#include "regionType.H"

#include "MULES.H"
#include "subCycle.H"
#include "UniformField.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{
    defineTypeNameAndDebug(driftFluxSystem, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftFluxSystem::driftFluxSystem
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),

    phase1_(phaseModels_[0]),
    phase2_(phaseModels_[1]),

    mixture_(new twoPhaseInteractingMixture(*this)),

    U_(mixture_->URef()),

    rho_
    (
        IOobject
        (
            "thermo:rho.mixture",
            mesh.time().timeName(),
            mesh
        ),
        mixture_->rho()
    ),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi.mixture",
            mesh.time().timeName(),
            mesh
        ),
        fvc::interpolate(rho_)*phi_
    ),

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
    ),

    UdmModel_(relativeVelocityModel::New(subDict("mixture"), mixture_()))
{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);

    phase2_.max(phase2_.residualAlpha().value());

    volScalarField& alpha1 = phase1_;
    mesh.schemes().setFluxRequired(alpha1.name());

    //- Update based on averaged velocity
    phi_ = fvc::flux(U_);

    rhoPhi_ = fvc::interpolate(rho_) * phi_;

    turbulence_ = compressible::momentumTransportModel
        ::New(rho_, U_, rhoPhi_, mixture_());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftFluxSystem::~driftFluxSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::driftFluxSystem::phirMag() const
{
    return mag(phi_);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::driftFluxSystem::phiMagMax() const
{
    return mag(phi_);
}


bool Foam::driftFluxSystem::isSinglePhase() const
{
    return false;
}


void Foam::driftFluxSystem::solve()
{
    Info << "\nSolving for mixture phase flow:" << endl;

    fvMesh& mesh = const_cast<fvMesh&>(mesh_);

    const Time& runTime = mesh.time();

    pimpleControl pimple(mesh);

    #include "createFields.H"

    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    #include "pUf/createRDeltaTf.H"

    if (LTS && faceMomentum)
    {
        #include "./setRDeltaTf.H"
    }

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        #include "alphaEqn.H"

        correct();

        mixture_->correct();

        #include "YEqns.H"

        if (faceMomentum)
        {
            #include "pUf/UEqns.H"
            #include "EEqns.H"
            #include "pUf/pEqn.H"
        }
        else
        {
            #include "pU/UEqns.H"
            #include "EEqns.H"
            #include "pU/pEqn.H"
        }

        correctKinematics();

        if (pimple.turbCorr())
        {
            turbulence_->correct();
        }
    }
}

// ************************************************************************* //
