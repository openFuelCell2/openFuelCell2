/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | openFuelCell
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

#include "addToRunTimeSelectionTable.H"
#include "powerLaw.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porousZones
    {
        defineTypeNameAndDebug(powerLaw, 0);
        addToRunTimeSelectionTable(porousZone, powerLaw, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZones::powerLaw::powerLaw
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porousZone(name, modelType, mesh, dict, cellZoneName),
    C0_(readScalar(coeffs_.lookup("C0"))),
    C1_(readScalar(coeffs_.lookup("C1"))),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porousZones::powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZones::powerLaw::calcTranformModelData()
{
    // nothing to be transformed
}


void Foam::porousZones::powerLaw::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    const scalarField& V = mesh_.V();

    apply(Udiag, V, rho, U);

    force = Udiag*U;
}


void Foam::porousZones::powerLaw::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
 
    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(Udiag, V, rho, U);
    }
    else
    {
        apply(Udiag, V, geometricOneField(), U);
    }
}


void Foam::porousZones::powerLaw::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
 
    apply(Udiag, V, rho, U);
}


void Foam::porousZones::powerLaw::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(AU, rho, U);
    }
    else
    {
        apply(AU, geometricOneField(), U);
    }
}

void Foam::porousZones::powerLaw::dpcds
(
    volTensorField& pc,
    volTensorField& dpcds
)
{
    //
}

void Foam::porousZones::powerLaw::correctU
(
    volVectorField& U1,
    volVectorField& U2,
    const volVectorField& U,
    const volTensorField& dpcdc
)
{
    // do nothing
}


bool Foam::porousZones::powerLaw::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
