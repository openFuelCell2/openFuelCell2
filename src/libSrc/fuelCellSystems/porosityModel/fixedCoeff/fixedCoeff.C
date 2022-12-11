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

#include "addToRunTimeSelectionTable.H"
#include "fixedCoeff.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porousZones
    {
        defineTypeNameAndDebug(fixedCoeff, 0);
        addToRunTimeSelectionTable(porousZone, fixedCoeff, mesh);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porousZones::fixedCoeff::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const vectorField& U,
    const scalar rho
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const label j = fieldIndex(i);
            const tensor Cd = rho*(alphaZones[j] + betaZones[j]*mag(U[cellI]));
            const scalar isoCd = tr(Cd);

            Udiag[cellI] += V[cellI]*isoCd;
            Usource[cellI] -= V[cellI]*((Cd - I*isoCd) & U[cellI]);
        }
    }
}


void Foam::porousZones::fixedCoeff::apply
(
    tensorField& AU,
    const vectorField& U,
    const scalar rho
) const
{

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const label j = fieldIndex(i);
            const tensor alpha = alphaZones[j];
            const tensor beta = betaZones[j];

            AU[cellI] += rho*(alpha + beta*mag(U[cellI]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZones::fixedCoeff::fixedCoeff
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porousZone(name, modelType, mesh, dict, cellZoneName),
    alphaXYZ_("alpha", dimless/dimTime, coeffs_),
    betaXYZ_("beta", dimless/dimLength, coeffs_),
    alpha_(cellZoneIDs_.size()),
    beta_(cellZoneIDs_.size())
{
    adjustNegativeResistance(alphaXYZ_);
    adjustNegativeResistance(betaXYZ_);

    calcTranformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porousZones::fixedCoeff::~fixedCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZones::fixedCoeff::calcTranformModelData()
{
    // The alpha coefficient as a tensor                   // S.Hess
    tensor alphaCoeff(Zero);
    alphaCoeff.xx() = alphaXYZ_.value().x();
    alphaCoeff.yy() = alphaXYZ_.value().y();
    alphaCoeff.zz() = alphaXYZ_.value().z();

    // The beta coefficient as a tensor
    tensor betaCoeff(Zero);
    betaCoeff.xx() = betaXYZ_.value().x();
    betaCoeff.yy() = betaXYZ_.value().y();
    betaCoeff.zz() = betaXYZ_.value().z();

    if (coordSys_.uniform())
    {
        forAll(cellZoneIDs_, zonei)
        {
            alpha_[zonei].resize(1);
            beta_[zonei].resize(1);

            alpha_[zonei] = coordSys_.transform(alphaCoeff);
            beta_[zonei] = coordSys_.transform(betaCoeff);
        }
    }
    else
    {
        forAll(cellZoneIDs_, zonei)
        {
            const pointUIndList cc
            (
                mesh_.cellCentres(),
                mesh_.cellZones()[cellZoneIDs_[zonei]]
            );

            alpha_[zonei] = coordSys_.transform(cc, alphaCoeff);
            beta_[zonei] = coordSys_.transform(cc, betaCoeff);
        }
    }
}


void Foam::porousZones::fixedCoeff::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), vector::zero);
    const scalarField& V = mesh_.V();
    scalar rhoRef = readScalar(coeffs_.lookup("rhoRef"));

    apply(Udiag, Usource, V, U, rhoRef);

    force = Udiag*U - Usource;
}


void Foam::porousZones::fixedCoeff::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porousZones::fixedCoeff::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porousZones::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(AU, U, rho);
}

void Foam::porousZones::fixedCoeff::dpcds
(
    volTensorField& pc,
    volTensorField& dpcds
)
{
    //
}


void Foam::porousZones::fixedCoeff::correctU
(
    volVectorField& U1,
    volVectorField& U2,
    const volVectorField& U,
    const volTensorField& dpcdc
)
{
    // do nothing
}

bool Foam::porousZones::fixedCoeff::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
