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
#include "electric.H"
#include "IOdictionary.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "sigmaModelList.H"
#include "dissolvedModel.H"
#include "activationOverpotentialModel.H"

#include "fuelCellSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(electric, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        electric,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::electric::electric
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    i_
    (
        IOobject
        (
            "i",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedVector(dimCurrent/dimArea, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    phi_
    (
        IOobject
        (
            activationOverpotentialModel::phiName,
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    j_
    (
        IOobject
        (
            activationOverpotentialModel::jName,
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimCurrent/dimVol, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    jPrime_
    (
        IOobject
        (
            "jPrime",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions()/phi_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    jSource_
    (
        IOobject
        (
            "jSource",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaField_
    (
        IOobject
        (
            "sigma",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions()/phi_.dimensions()*dimArea, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimTemperature, dict_.getOrDefault<scalar>("T0", 293.15)),
        zeroGradientFvPatchScalarField::typeName
    ),
    control_(dict_.lookupOrDefault<Switch>("control", false)),
    dissolveOnOff_(dict_.lookupOrDefault<Switch>("dissolveOnOff", false)),
    patchName_(word::null),
    zoneName_(word::null),
    cellZoneIDs_(),
    galvanostatic_(true),
    ibar_(nullptr),
    voltage_(nullptr),
    maxPhiStep_(0.2)
{
    if (dissolveOnOff_)
    {
        dissolved_ = dissolvedModel::New(*this, dict_.subDict("dissolved"));
    }

    sigma_.reset
    (
        new sigmaModelList
        (
            *this,
            dict_.subDict("sigma")
        )
    );

    if (control_)
    {
        const dictionary& controlDict = dict_.subDict("galvanostatic");
        patchName_ = controlDict.get<word>("patchName");
        galvanostatic_ = controlDict.get<Switch>("active");

        if (galvanostatic_)
        {
            ibar_.reset(Function1<scalar>::New("ibar", controlDict, this));
            maxPhiStep_ =
                controlDict.lookupOrDefault<scalar>("maxPhiStep", 0.2);

            if (maxPhiStep_ <= 0)
            {
                FatalIOErrorInFunction(controlDict)
                    << "maxPhiStep must be greater than zero"
                    << exit(FatalIOError);
            }
        }
        else
        {
            voltage_.reset(Function1<scalar>::New("voltage", controlDict, this));
        }

        if (controlDict.found("zoneName"))
        {
            zoneName_ = controlDict.get<word>("zoneName");
            cellZoneIDs_ = this->cellZones().indices(zoneName_);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::electric::~electric()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::electric::solvePotential()
{
    tmp<fvScalarMatrix> phiEqn
    (
      - fvm::laplacian(sigmaField_, phi_, "laplacian(sigma,phi)")
      + fvm::Sp(jPrime_, phi_)
     == jSource_
    );

    // A single reference is sufficient.  The positive Newton diagonal
    // normally removes the Laplacian null space in reacting regions.
    if (phi_.needReference() && nCells())
    {
        phiEqn->setReference(0, phi_[0], true);
    }

    phiEqn->solve();
}


Foam::scalar Foam::regionTypes::electric::integratedCurrent() const
{
    const scalarField source
    (
        jSource_.primitiveField()
      - jPrime_.primitiveField()*phi_.primitiveField()
    );

    scalar current = 0;

    if (cellZoneIDs_.size())
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelList& cells = cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                current += source[cellI]*V()[cellI];
            }
        }
        reduce(current, sumOp<scalar>());
    }
    else
    {
        current = Foam::gSum(source*V());
    }

    return current;
}


void Foam::regionTypes::electric::solve()
{
    Info << "\nSolve for region " << name() << ":\n" << endl;

    if (control_ && galvanostatic_)
    {
        const label patchID = boundaryMesh().findPatchID(patchName_);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Cannot find " << patchName_
                << "; please check patchName" << exit(FatalError);
        }

        const scalar area = Foam::gSum(magSf().boundaryField()[patchID]);
        const scalar oldPhi = Foam::gAverage(phi_.boundaryField()[patchID]);

        // The Newton-linearized equation is affine in the controlled
        // boundary value.  Two unit-response solves therefore give the
        // boundary potential required for the prescribed total current.
        phi_.boundaryFieldRef()[patchID] == scalar(0);
        solvePotential();
        const scalar current0 = integratedCurrent();
        const scalarField phi0(phi_.primitiveField());

        phi_.boundaryFieldRef()[patchID] == scalar(1);
        solvePotential();
        const scalar current1 = integratedCurrent();
        const scalarField phi1(phi_.primitiveField());

        const scalar conductance = current1 - current0;
        const scalar targetCurrent = ibar_->value(time().value())*area;
        scalar targetPhi = oldPhi;

        if (mag(conductance) > VSMALL)
        {
            targetPhi = (targetCurrent - current0)/conductance;
        }
        else
        {
            WarningInFunction
                << "The galvanostatic unit response is too small in region "
                << name() << "; retaining boundary potential " << oldPhi
                << endl;
        }

        const scalar deltaPhi = Foam::max
        (
            -maxPhiStep_,
            Foam::min(maxPhiStep_, targetPhi - oldPhi)
        );
        const scalar boundaryPhi = oldPhi + deltaPhi;

        phi_.primitiveFieldRef() =
            phi0 + boundaryPhi*(phi1 - phi0);
        phi_.boundaryFieldRef()[patchID] == boundaryPhi;
        phi_.correctBoundaryConditions();

        const scalar current = current0 + boundaryPhi*conductance;
        Info<< "ibar: " << current/area
            << '\t' << "voltage: " << boundaryPhi << endl;
    }
    else
    {
        solvePotential();

        if (control_)
        {
            const label patchID = boundaryMesh().findPatchID(patchName_);
            const scalar area = Foam::gSum(magSf().boundaryField()[patchID]);
            Info<< "ibar: " << integratedCurrent()/area
                << '\t' << "voltage: "
                << Foam::gAverage(phi_.boundaryField()[patchID]) << endl;
        }
    }

    i_ = -sigmaField_ * fvc::grad(phi_);
    i_.correctBoundaryConditions();

    if (dissolved_.valid())
    {
        dissolved_->solve();
    }
}


void Foam::regionTypes::electric::setRDeltaT()
{
    //- Do nothing, add if necessary.
}


void Foam::regionTypes::electric::correct()
{
    sigma_->correct(sigmaField_);
    sigmaField_.correctBoundaryConditions();

    // Models that do not provide a Jacobian retain the original explicit
    // source treatment.
    if (Foam::gMax(jPrime_.primitiveField()) <= VSMALL)
    {
        jSource_ = j_;
    }

    if (control_ && !galvanostatic_)
    {
        const label patchID = boundaryMesh().findPatchID(patchName_);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Cannot find " << patchName_
                << "; please check patchName" << exit(FatalError);
        }

        phi_.boundaryFieldRef()[patchID] == voltage_->value(time().value());
    }

    if (dissolved_.valid())
    {
        dissolved_->correct();
    }
}


void Foam::regionTypes::electric::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " to Cell " << nl << endl;

    //- heat source
    volScalarField heatSource = (i_ & i_)/sigmaField_;
    volScalarField heatSource0
    (
        IOobject
        (
            "heatSource",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    heatSource0.rmap(heatSource, cellMapIO_);

    fuelCell.Qdot() += heatSource0;
}


void Foam::regionTypes::electric::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " from Cell " << nl << endl;

    scalarField& T = T_;

    forAll(T, cellI)
    {
        T[cellI] = fuelCell.T()[cellMapIO_[cellI]];
    }

    T_.correctBoundaryConditions();
}

// ************************************************************************* //
