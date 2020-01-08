/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
#include "electronIon.H"
#include "IOdictionary.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "sigmaModelList.H"
#include "dissolvedModel.H"

#include "fuelCellSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(electronIon, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        electronIon,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::electronIon::electronIon
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regionType(mesh, dict),

    mesh_(mesh),
    dict_(dict),
    i_
    (
        IOobject
        (
            "i",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedVector("i", dimensionSet(0, -2, 0, 0, 0, 1, 0), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    phi_
    (
        IOobject
        (
            "phi",
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
            "J",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("J", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaField_
    (
        IOobject
        (
            "sigma",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("sigmaE", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("T", dimTemperature, 353.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    relax_(dict_.lookupOrDefault<scalar>("relax", 0.0)),
    active_(false),
    patchName_(word::null),
    galvanostatic_(true),
    ibar_(0.0),
    voltage_(0.0)
{
    const dictionary& cellProperties =
        mesh_.lookupObject<IOdictionary>("cellProperties");

    const dictionary& dissolvedDict = cellProperties.subDict("dissolved");

    if
    (
        word(dissolvedDict.lookup("region"))
     == this->name()
    )
    {
        dissolved_ = dissolvedModel::New(*this, dissolvedDict);
    }

    sigma_.set
    (
        new sigmaModelList
        (
            *this,
            dict_.subDict("sigma")
        )
    );

    const dictionary& galvanostaticDict = cellProperties.subDict("galvanostatic");

    if
    (
        word(galvanostaticDict.lookup("region"))
     == this->name()
    )
    {
        active_ = true;
        galvanostatic_ = Switch(galvanostaticDict.lookupOrDefault("active", true));
        patchName_ = word(galvanostaticDict.lookup("patchName"));

        if (galvanostatic_)
        {
            galvanostaticDict.lookup("ibar") >> ibar_;
        }
        else
        {
            galvanostaticDict.lookup("voltage") >> voltage_;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::electronIon::~electronIon()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::electronIon::solve()
{
    tmp<fvScalarMatrix> phiEqn
    (
      - fvm::laplacian(sigmaField_, phi_)
      - j_
    );

    phiEqn->relax();

    phiEqn->solve();

    i_ = -sigmaField_ * fvc::grad(phi_);
    i_.correctBoundaryConditions();

    if (dissolved_.valid())
    {
        dissolved_->solve();
    }
}


void Foam::regionTypes::electronIon::setRDeltaT()
{
    //- Do nothing, add if necessary.
}


void Foam::regionTypes::electronIon::correct()
{
    sigma_->correct(sigmaField_);
    sigmaField_.correctBoundaryConditions();

    if (phi_.needReference() || active_)
    {
        //- Update the potential field.
        //- ElectroNatural, the total electron+proton flux should be zero
        const scalarField& source = j_;
        const scalarField& volume = this->V();
        scalarField sum = source*volume;
        scalar iDot(Foam::gSum(sum));

        if (active_)
        {
            label patchID = this->boundaryMesh().findPatchID(patchName_);

            if (patchID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find " << patchName_
                    << "please check the patchName" << exit(FatalError);
            }

            scalarField phiBoundary = phi_.boundaryFieldRef()[patchID];

            scalar ibar0 = iDot/Foam::gSum(this->magSf().boundaryField()[patchID]);

            phi_.boundaryFieldRef()[patchID] == phiBoundary + relax_*(ibar0 - ibar_);

            Info << "ibar: " << ibar0 << endl;
            Info << "voltage: " << phiBoundary[0] << endl;
        }
        else
        {
            scalarField& phi = phi_;
            phi -= iDot*relax_;
        }
    }

    if (dissolved_.valid())
    {
        dissolved_->correct();
    }
}


void Foam::regionTypes::electronIon::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map" << name() << " to Cell" << endl;

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
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );

    heatSource0.rmap(heatSource, cellMapIO_);

    fuelCell.Qdot() += heatSource0;
}


void Foam::regionTypes::electronIon::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " from Cell" << endl;

    scalarField& T = T_;

    forAll(T, cellI)
    {
        T[cellI] = fuelCell.T()[cellMapIO_[cellI]];
    }

    T_.correctBoundaryConditions();
}

// ************************************************************************* //
