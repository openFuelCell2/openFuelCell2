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
#include "fuelCellSystem.H"
#include "phaseSystem.H"
#include "HashPtrTable.H"
#include "patchToPatchInterpolation.H"
#include "regionTypeList.H"


// * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * //

Foam::IOobject Foam::fuelCellSystem::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "cellProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Get cell properties from " << io.name() << nl << endl;

        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        Info<< "No cell properties are provided... " << nl << endl;

        io.readOpt(IOobject::NO_READ);
    }

    return io;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fuelCellSystem::fuelCellSystem
(
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),

    mesh_(mesh),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Qdot_
    (
        T_,
        dimEnergy/dimTime
    ),

    k_
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    contErrCp_
    (
        IOobject
        (
            "contErrCp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimDensity/dimTime*dimSpecificHeatCapacity,
            Zero
        )
    ),

    Cp_
    (
        IOobject
        (
            "Cp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimSpecificHeatCapacity,
            Zero
        )
    ),

    rhoCp_
    (
        IOobject
        (
            "rhoCp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimDensity*dimSpecificHeatCapacity,
            Zero
        )
    ),

    rhoCpPhi_
    (
        IOobject
        (
            "rhoCpPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimVelocity*dimDensity*dimSpecificHeatCapacity*dimArea,
            Zero
        )
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVelocity*dimArea, Zero)
    )
{
    Info << nl << "Creating regions: " << endl;
    regions_.set
    (
        new regionTypeList
        (
            mesh_
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fuelCellSystem::~fuelCellSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fuelCellSystem::correct()
{
    //- regions
    regions_->correct();
}


void Foam::fuelCellSystem::setRDeltaT()
{
    regions_->setRDeltaT();
}


void Foam::fuelCellSystem::solve()
{
    //- solve for regions
    regions_->solve();
}


void Foam::fuelCellSystem::mapToCell()
{
    //- T source, = 0

    Qdot_ *= 0.0;

    Info<< "\nMap to cell \n" << endl;

    regions_->mapToCell(*this);
}


void Foam::fuelCellSystem::mapFromCell()
{
    Info<< "\nMap from cell \n" << endl;

    regions_->mapFromCell(*this);
}


Foam::regionTypeList& Foam::fuelCellSystem::regions()
{
    return regions_();
}

// ************************************************************************* //
