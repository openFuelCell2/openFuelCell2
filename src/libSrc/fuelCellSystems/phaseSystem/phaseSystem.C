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
#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "aspectRatioModel.H"
#include "surfaceInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "dragModel.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
    defineRunTimeSelectionTable(phaseSystem, dictionary);
}

const Foam::word Foam::phaseSystem::propertiesName("regionProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
Foam::typeIOobject<Foam::IOdictionary>
Foam::phaseSystem::readPhasePropertiesDict
(
    const objectRegistry& obr
)
{
    typeIOobject<IOdictionary> phasePropertiesIO
    (
        propertiesName,
        obr.time().constant(),
        obr,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        true
    );

    if (phasePropertiesIO.headerOk())
    {
        Info<< "Get phase properties from " << phasePropertiesIO.name() << nl << endl;
    }
    else
    {
        Info<< "No phase properties presented..." << nl << endl;

        phasePropertiesIO.readOpt() = IOobject::NO_READ;
    }

    return phasePropertiesIO;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::calcPhi
(
    const phaseModelList& phaseModels
) const
{
    tmp<surfaceScalarField> tmpPhi
    (
        new surfaceScalarField
        (
            "phi",
            fvc::interpolate(phaseModels[0])*phaseModels[0].phi()
        )
    );

    for (label phasei=1; phasei<phaseModels.size(); phasei++)
    {
        tmpPhi.ref() +=
            fvc::interpolate(phaseModels[phasei])*phaseModels[phasei].phi();
    }

    return tmpPhi;
}


void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {}

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary(readPhasePropertiesDict(mesh)),

    mesh_(mesh),

    porousZone_
    (
        new porousZoneList
        (
            mesh
        )
    ),

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    phi_(calcPhi(phaseModels_)),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", dimPressure/dimTime, 0)
    ),

    MRF_(mesh_),

    continuous_(lookupOrDefault<word>("continuous", phaseModels_[0].name())),

    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        dimensionedVector(dimVelocity/dimTime, vector::zero)
    ),

    pRefCell_(0),
    pRefValue_(0.0),
    pimple_(const_cast<fvMesh&>(mesh))
{
    // Pimple control
    setRefCell
    (
        phaseModels_[0].thermoRef().p(),
        p_rgh_,
        pimple_.dict(),
        pRefCell_,
        pRefValue_
    );
    mesh.schemes().setFluxRequired(p_rgh_.name());

    // Create LTS field, if LTS is enabled
    bool LTS = fv::localEulerDdt::enabled(mesh);

    if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT_ = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("one", dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }

    //- Get gravity field
    {
        uniformDimensionedVectorField g
        (
            IOobject
            (
                "g",
                mesh.time().constant(),
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        g_[0] = g[0];
    }

    // Groupings
    label movingPhasei = 0;
    label stationaryPhasei = 0;
    label anisothermalPhasei = 0;
    label multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        movingPhasei += !phase.stationary();
        stationaryPhasei += phase.stationary();
        anisothermalPhasei += !phase.isothermal();
        multiComponentPhasei += !phase.pure();
    }
    movingPhaseModels_.resize(movingPhasei);
    stationaryPhaseModels_.resize(stationaryPhasei);
    anisothermalPhaseModels_.resize(anisothermalPhasei);
    multiComponentPhaseModels_.resize(multiComponentPhasei);

    movingPhasei = 0;
    stationaryPhasei = 0;
    anisothermalPhasei = 0;
    multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        if (!phase.stationary())
        {
            movingPhaseModels_.set(movingPhasei ++, &phase);
        }
        if (phase.stationary())
        {
            stationaryPhaseModels_.set(stationaryPhasei ++, &phase);
        }
        if (!phase.isothermal())
        {
            anisothermalPhaseModels_.set(anisothermalPhasei ++, &phase);
        }
        if (!phase.pure())
        {
            multiComponentPhaseModels_.set(multiComponentPhasei ++, &phase);
        }
    }

    // Creation only for multiphase flow.
    if(phaseModels_.size() > 1)
    {
        // Blending methods
        forAllConstIter(dictionary, subDict("blending"), iter)
        {
            blendingMethods_.insert
            (
                iter().keyword(),
                blendingMethod::New
                (
                    iter().keyword(),
                    iter().dict(),
                    phaseModels_.toc()
                )
            );
        }

        // Sub-models
        generatePairsAndSubModels("surfaceTension", surfaceTensionModels_);
        generatePairsAndSubModels("aspectRatio", aspectRatioModels_);
    }

    // Write phi
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    // Update motion fields
    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    const label nMovingPhases = movingPhaseModels_.size();

    tmp<volScalarField> rho(movingPhaseModels_[0]*movingPhaseModels_[0].rho());
    for (label movingPhasei = 1; movingPhasei < nMovingPhases; ++ movingPhasei)
    {
        rho.ref() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].rho();
    }

    if (stationaryPhaseModels_.empty())
    {
        return rho;
    }

    volScalarField alpha(movingPhaseModels_[0]);
    for (label movingPhasei = 1; movingPhasei < nMovingPhases; ++ movingPhasei)
    {
        alpha += movingPhaseModels_[movingPhasei];
    }

    return rho/alpha;
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    const label nMovingPhases = movingPhaseModels_.size();

    tmp<volVectorField> U(movingPhaseModels_[0]*movingPhaseModels_[0].U());
    for (label movingPhasei = 1; movingPhasei < nMovingPhases; ++ movingPhasei)
    {
        U.ref() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].U();
    }

    if (stationaryPhaseModels_.empty())
    {
        return U;
    }

    volScalarField alpha(movingPhaseModels_[0]);
    for (label movingPhasei = 1; movingPhasei < nMovingPhases; ++ movingPhasei)
    {
        alpha += movingPhaseModels_[movingPhasei];
    }

    return U/alpha;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    aspectRatioModel::typeName + ":E",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", dimless, 1)
            )
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phasePairKey& key) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma();
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    surfaceTensionModel::typeName + ":sigma",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", surfaceTensionModel::dimSigma, 0)
            )
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::dmdt
(
    const phasePairKey& key
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dmdt", phasePairs_[key]->name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );
}


Foam::PtrList<Foam::volScalarField> Foam::phaseSystem::dmdts() const
{
    PtrList<volScalarField> dmdts(this->phaseModels_.size());

    return dmdts;
}


void Foam::phaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}


void Foam::phaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctKinematics();

        updateDpdt = updateDpdt || phaseModels_[phasei].thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.begin()().thermo().p());
    }
}


void Foam::phaseSystem::correctThermo()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermo();
    }
}


void Foam::phaseSystem::correctTurbulence()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctTurbulence();
    }
}


void Foam::phaseSystem::correctEnergyTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctEnergyTransport();
    }
}

void Foam::phaseSystem::correctBoundaryFlux()
{
    forAll(movingPhases(), movingPhasei)
    {
        phaseModel& phase = movingPhases()[movingPhasei];

        const volVectorField::Boundary& UBf = phase.U()().boundaryField();

        FieldField<fvsPatchField, scalar> phiRelBf
        (
            MRF_.relative(mesh_.Sf().boundaryField() & UBf)
        );

        surfaceScalarField::Boundary& phiBf = phase.phiRef().boundaryFieldRef();

        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                isA<fixedValueFvsPatchScalarField>(phiBf[patchi])
             && !isA<movingWallVelocityFvPatchVectorField>(UBf[patchi])
            )
            {
                phiBf[patchi] == phiRelBf[patchi];
            }
        }
    }
}


Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::phaseSystem::heatTransfer() const
{
    autoPtr<heatTransferTable> eqnsPtr
    (
        new heatTransferTable()
    );

    heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.set
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    return eqnsPtr;
}


Foam::autoPtr<Foam::fvScalarMatrix>
Foam::phaseSystem::heatTransfer
(
    volScalarField& T,
    labelList& cellMap
) const
{
    return autoPtr<fvScalarMatrix>(new fvScalarMatrix(T, dimEnergy/dimTime));
}


void Foam::phaseSystem::setRDeltaT()
{
    if (!trDeltaT_.valid())
    {
        return;
    }

    volScalarField& rDeltaT = trDeltaT_();

    const dictionary& LTSdict = this->subDict("LTS");

    scalar maxCo
    (
        LTSdict.lookupOrDefault<scalar>("maxCo", 0.2)
    );

    scalar maxDeltaT
    (
        LTSdict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    scalar rDeltaTSmoothingCoeff
    (
        LTSdict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(phiMagMax())()()
       /((2*maxCo)*mesh_.V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        forAll(phaseModels_, phasei)
        {
            readOK &= phaseModels_[phasei].read();
        }

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField> Foam::byDt(const volScalarField& vf)
{
    if (fv::localEulerDdt::enabled(vf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaT(vf.mesh())*vf;
    }
    else
    {
        return vf/vf.mesh().time().deltaT();
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::byDt(const surfaceScalarField& sf)
{
    if (fv::localEulerDdt::enabled(sf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaTf(sf.mesh())*sf;
    }
    else
    {
        return sf/sf.mesh().time().deltaT();
    }
}


// ************************************************************************* //
