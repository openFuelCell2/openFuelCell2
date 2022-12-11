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
#include "addToRunTimeSelectionTable.H"
#include "DarcyForchheimer.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porousZones
    {
        defineTypeNameAndDebug(DarcyForchheimer, 0);
        addToRunTimeSelectionTable(porousZone, DarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZones::DarcyForchheimer::DarcyForchheimer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porousZone(name, modelType, mesh, dict, cellZoneName),
    dXYZ_("d", dimless/sqr(dimLength), coeffs_),
    fXYZ_("f", dimless/dimLength, coeffs_),
    D_(cellZoneIDs_.size()),
    sqrtD_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "thermo:rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "thermo:nu")),
    alphaName_(coeffs_.lookupOrDefault<word>("alpha", "alpha")),
    theta_(coeffs_.lookupOrDefault<scalar>("theta", 0.0)),
    sigma_(coeffs_.lookupOrDefault<scalar>("sigma", 0.0))
{
    adjustNegativeResistance(dXYZ_);
    adjustNegativeResistance(fXYZ_);

    calcTranformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porousZones::DarcyForchheimer::~DarcyForchheimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZones::DarcyForchheimer::calcTranformModelData()
{
    // The Darcy coefficient as a tensor           // S.Hess
    tensor darcyCoeff(Zero);
    darcyCoeff.xx() = dXYZ_.value().x();
    darcyCoeff.yy() = dXYZ_.value().y();
    darcyCoeff.zz() = dXYZ_.value().z();

    // The Forchheimer coefficient as a tensor
    // - the leading 0.5 is from 1/2*rho
    tensor forchCoeff(Zero);
    forchCoeff.xx() = 0.5*fXYZ_.value().x();
    forchCoeff.yy() = 0.5*fXYZ_.value().y();
    forchCoeff.zz() = 0.5*fXYZ_.value().z();

    tensor sqrtDCoeff(Zero);
    sqrtDCoeff.xx() = sqrt(dXYZ_.value().x());
    sqrtDCoeff.yy() = sqrt(dXYZ_.value().y());
    sqrtDCoeff.zz() = sqrt(dXYZ_.value().z());

    if (coordSys_.uniform())
    {
        forAll(cellZoneIDs_, zonei)
        {
            D_[zonei].resize(1);
            F_[zonei].resize(1);
            sqrtD_[zonei].resize(1);                // S.Hess

            sqrtD_[zonei] = coordSys_.transform(sqrtDCoeff);    // S.Hess
            D_[zonei] = coordSys_.transform(darcyCoeff);
            F_[zonei] = coordSys_.transform(forchCoeff);
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

            D_[zonei] = coordSys_.transform(cc, darcyCoeff);
             sqrtD_[zonei] = coordSys_.transform(cc, sqrtDCoeff);
            F_[zonei] = coordSys_.transform(cc, forchCoeff);
        }
    }                                                           // up to here

    if (debug && mesh_.time().outputTime())
    {
        volTensorField Dout
        (
            IOobject
            (
                typeName + ":D",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("0", dXYZ_.dimensions(), tensor::zero)
        );
        volTensorField Fout
        (
            IOobject
            (
                typeName + ":F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("0", fXYZ_.dimensions(), tensor::zero)
        );

        UIndirectList<tensor>(Dout, mesh_.cellZones()[cellZoneIDs_[0]]) = D_[0];
        UIndirectList<tensor>(Fout, mesh_.cellZones()[cellZoneIDs_[0]]) = F_[0];

        Dout.write();
        Fout.write();
    }
}


void Foam::porousZones::DarcyForchheimer::calcForce
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

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;
}


void Foam::porousZones::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));
    word alphaName(IOobject::groupName(alphaName_, U.group()));

    //- alpha field
    tmp<volScalarField> alpha
    (
        new volScalarField
        (
            IOobject
            (
                "alpha0",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimless, 1.0)
        )
    );

    if (mesh_.foundObject<volScalarField>(alphaName))
    {
        alpha.ref() = mesh_.lookupObject<volScalarField>(alphaName);
    }
    else
    {
        //- look up the full name of alpha field
        word alphaFullName(coeffs_.lookup("alphaFull"));

        alpha.ref() = mesh_.lookupObject<volScalarField>(alphaFullName);
    }

    alpha.ref().max(Foam::SMALL); //avoid

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu/alpha, U);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu/alpha, U);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu/alpha, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho/alpha, U);
        }
    }
}


void Foam::porousZones::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, mu, U);
}


void Foam::porousZones::DarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU, rho, mu, U);
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(AU, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(AU, geometricOneField(), mu/rho, U);
        }
    }
}


void Foam::porousZones::DarcyForchheimer::dpcds
(
    volTensorField& pc,
    volTensorField& dpcds
)
{
    word phaseName(coeffs_.lookupOrDefault<word>("water", "water"));

    word alphaName(IOobject::groupName(alphaName_, phaseName));

    if(mesh_.foundObject<volScalarField>(alphaName))
    {
        const volScalarField& alpha = mesh_.lookupObject<volScalarField>(alphaName);

        forAll(cellZoneIDs_, zoneI)
        {
            const tensorField& dZones = sqrtD_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                const label j = this->fieldIndex(i);
                pc[cellI]=
                    sigma_*cos(theta_)*sqrt(porosity_)*dZones[j]
                    *(
                        theta_ < Foam::constant::mathematical::pi/2.0
                        ? (1.42*(1.0 - alpha[cellI]) - 2.12*pow(1.0 - alpha[cellI], 2) + 1.26*pow(1.0 - alpha[cellI], 3))
                        : (1.42*alpha[cellI] - 2.12*pow(alpha[cellI], 2) + 1.26*pow(alpha[cellI], 3))
                    );

                dpcds[cellI] =
                    sigma_*cos(theta_)*sqrt(porosity_)*dZones[j]
                    *(
                        theta_ < Foam::constant::mathematical::pi/2.0
                        ? (-1.42 + 2.0*2.12*(1.0 - alpha[cellI]) - 3.0*1.26*pow(1.0 - alpha[cellI], 2))
                        : (1.42 -2.0*2.12*alpha[cellI] + 3.0*1.26*pow(alpha[cellI], 2))
                    );
            }
        }
    }
}


void Foam::porousZones::DarcyForchheimer::correctU
(
    volVectorField& U1,
    volVectorField& U2,
    const volVectorField& U,
    const volTensorField& dpcds
)
{
    word mu1Name(IOobject::groupName(muName_, U1.group()));
    word alpha1Name(IOobject::groupName(alphaName_, U1.group()));
    word mu2Name(IOobject::groupName(muName_, U2.group()));
    word alpha2Name(IOobject::groupName(alphaName_, U2.group()));

    //- Dynamic viscosity
    const volScalarField& mu1 = mesh_.lookupObject<volScalarField>(mu1Name);
    const volScalarField& mu2 = mesh_.lookupObject<volScalarField>(mu2Name);

    //- phase saturation
    const volScalarField& alpha1 = mesh_.lookupObject<volScalarField>(alpha1Name);
    const volScalarField& alpha2 = mesh_.lookupObject<volScalarField>(alpha2Name);

    //- temperary velocities
    volVectorField U10 = dpcds & fvc::grad(alpha1);

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const label j = this->fieldIndex(i);

            //- only update velocity in porous media
            U2[cellI] = U[cellI];

            U1[cellI] = - Foam::pow(alpha1[cellI], 2)/mu1[cellI]*U10[cellI]/dZones[j]
                        + (
                              Foam::pow(alpha2[cellI], 2) * mu2[cellI]
                            / Foam::pow(alpha1[cellI], 2) / mu1[cellI]
                          ) * U[cellI];  
        }
    }
}



bool Foam::porousZones::DarcyForchheimer::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
