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
#include "fluid.H"
#include "fuelCellSystem.H"
#include "patchToPatchInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(fluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        fluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::fluid::fluid
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh)
{
    //- correct fluid, solve the mass and momentum equations, singlePhase and twoPhase.
    phases_ = phaseSystem::New
    (
        *this
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::regionTypes::fluid::~fluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::regionTypes::fluid::solve()
{
    phases_->solve();
}


void Foam::regionTypes::fluid::setRDeltaT()
{
    phases_->setRDeltaT();
}


void Foam::regionTypes::fluid::correct()
{
    phases_->correctEnergyTransport();
    phases_->correctThermo();
    phases_->correct();
}


void Foam::regionTypes::fluid::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " to Cell" << nl << endl;

    const uniformDimensionedVectorField& g =
        this->time().lookupObject<uniformDimensionedVectorField>("g");

    //- Continuous phase name
    const word& continuous = phases_->continuous();

    // Phase model
    phaseModel& phase = phases_->phases()[continuous];

    if (phase.isothermal())
    {
        return;
    }

    //- Gravity effect heat, alpha*rho*(U&g)
    volScalarField heatSource
    (
        phase
      * phase.thermo().rho()
      * (phase.U()&g)
    );

    forAll(phases_->phases(), phasei)
    {
        heatSource += phases_->phases()[phasei].Qdot().ref();
    }

    // Heat transfer field in parent mesh
    volScalarField heatSource0
    (
        IOobject
        (
            "heatSource",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );

    heatSource0.rmap(heatSource, cellMapIO_);

    // Map to parent mesh
    fuelCell.Qdot() += heatSource0;
    fuelCell.Qdot() += phases_->heatTransfer(fuelCell.T(), cellMapIO_)();

    // Rho*Cp
    scalarField rhoCp
    (
        phase
      * phase.thermo().Cp()
      * phase.thermo().rho()
    );

    // contErr*Cp
    scalarField contErrCp
    (
        phase.continuityError()
      * phase.thermo().Cp()
    );

    // Perform reverse mapping
    fuelCell.rhoCp().rmap(rhoCp, cellMapIO_);

    fuelCell.contErrCp().rmap(contErrCp, cellMapIO_);

    // Map air fluxes
    labelList internalFaceMap
    (
        SubList<label>(faceMap_, this->nInternalFaces())
    );

    scalarField internalFaceMask
    (
        scalarField::subField(faceMask_, this->nInternalFaces())
    );

    //
    // ** recall phi already incorporates rho **
    //
    surfaceScalarField& rhoPhi = phase.alphaRhoPhiRef();

    scalarField rhoCpPhi
    (
        linearInterpolate(phase.thermo().Cp())
      * rhoPhi
    );

    fuelCell.phi().rmap
    (
        rhoPhi*internalFaceMask,
        internalFaceMap
    );

    fuelCell.rhoCpPhi().rmap
    (
        rhoCpPhi*internalFaceMask,
        internalFaceMap
    );

    // Do flux boundary conditions
    forAll (patchesMapIO_, patchI)
    {
        if
        (
            patchesMapIO_[patchI] > -1
         && patchesMapIO_[patchI] < mesh_.boundary().size()
        )
        {
            // Patch maps
            labelField curFpm
            (
                labelField::subField
                (
                    faceMap_,
                    this->boundary()[patchI].size(),
                    this->boundary()[patchI].patch().start()
                )
            );

            scalarField curMask
            (
                scalarField::subField
                (
                    faceMask_,
                    this->boundary()[patchI].size(),
                    this->boundary()[patchI].patch().start()
                )
            );

            curFpm -= mesh_.boundary()
                        [patchesMapIO_[patchI]].patch().start();

            fuelCell.phi().boundaryFieldRef()[patchesMapIO_[patchI]].
                scalarField::rmap
                (
                    (
                        rhoPhi.boundaryFieldRef()[patchI]
                    )*curMask,
                    curFpm
                );

            fuelCell.rhoCpPhi().boundaryFieldRef()[patchesMapIO_[patchI]].
                scalarField::rmap
                (
                    (
                        phase.thermo().Cp().ref().boundaryFieldRef()[patchI]
                      * rhoPhi.boundaryFieldRef()[patchI]
                    )*curMask,
                    curFpm
                );
        }
    }

    //- effective thermal conductivities

    scalarField kF(nCells(), 0.0);

    kF = phase.kappa()*phase;

    forAll(phases_->porousZone(), iz)
    {
        label znId =
                this->cellZones().findZoneID(phases_->porousZone()[iz].zoneName());

        scalar CpZn, kZn;
        phases_->porousZone()[iz].dict().lookup("Cp") >> CpZn;
        phases_->porousZone()[iz].dict().lookup("k") >> kZn;

        scalar porZn = phases_->porousZone()[iz].porosity();

        labelList znCells(this->cellZones()[znId]);

        forAll(znCells, cellI)
        {
            label cellId = znCells[cellI];

            kF[cellId] =
                    kZn*(scalar(1) - porZn) + kF[cellId]*porZn;
        }
    }

    // Perform reverse mapping
    fuelCell.k().rmap(kF, cellMapIO_);
}


void Foam::regionTypes::fluid::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " from Cell" << nl << endl;

    const word& continuous = phases_->continuous();

    phaseModel& phase = phases_->phases()[continuous];

    volScalarField& he = phase.thermoRef().he();
    volScalarField& p = phase.thermoRef().p();

    volScalarField T0 = phase.thermoRef().T();

    forAll(phase.thermo().T(), cellI)
    {
        T0[cellI] = fuelCell.T()[cellMapIO_[cellI]];
    }

    he = phase.thermoRef().he(p, T0).ref();

    he.correctBoundaryConditions();
}

// ************************************************************************* //
