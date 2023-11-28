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
#include "solid.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(solid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        solid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::solid::solid
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    Qdot_
    (
        IOobject
        (
            "Qdot",
            this->time().timeName(),
            *this
        ),
        *this,
        dimensionedScalar("Q", dimEnergy/dimTime/dimVolume, 0)
    )
{
    //- Solid thermo model
    thermo_ = solidThermo::New(*this);

    radiation_ = radiationModel::New(thermo_->T());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::solid::~solid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::solid::correct()
{
    thermo_->correct();
}


void Foam::regionTypes::solid::setRDeltaT()
{
    // do nothing, add something if necessary
}


void Foam::regionTypes::solid::solve()
{
    // do nothing, add something if necessary
}


void Foam::regionTypes::solid::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " to Cell" << nl << endl;

    volScalarField& T = thermo_->T();

    //- heat source
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

    const volScalarField rho = thermo_->rho();
    const volScalarField cp = thermo_->Cp();

    //- rhoCpCell
    volScalarField& rhoCpCell = fuelCell.rhoCp();

    scalarField rhoCp = rho*cp;

    scalarField radiationST =
        radiation_->Ru()/rhoCp
      - radiation_->Rp()*Foam::pow(T, 4)/rhoCp;

    heatSource0.rmap(Qdot_ + radiationST, cellMapIO_);

    fuelCell.Qdot() += heatSource0;

    // Perform reverse mapping
    rhoCpCell.rmap(rhoCp, cellMapIO_);

    rhoCpCell.correctBoundaryConditions();

    //- thermal conductivity
    volScalarField& kCell = fuelCell.k();
    const scalarField kappa = thermo_->kappa();

    kCell.rmap(kappa, cellMapIO_);
    kCell.correctBoundaryConditions();
}


void Foam::regionTypes::solid::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " from Cell" << nl << endl;

    volScalarField& he = thermo_->he();
    uniformGeometricScalarField& p = he.mesh().
        objectRegistry::lookupObjectRef<uniformGeometricScalarField>("p");

    // Define a pressure field
    tmp<volScalarField> p0 = volScalarField::New
    (
        "p0",
        *this,
        p
    );

    volScalarField T = thermo_->T();

    forAll(T, cellI)
    {
        T[cellI] = fuelCell.T()[cellMapIO_[cellI]];
    }

    he = thermo_->he(p0.ref(), T).ref();
}

// ************************************************************************* //
