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

#include "generalSigma.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{
    defineTypeNameAndDebug(generalSigma, 0);

    addToRunTimeSelectionTable
    (
        sigmaModel,
        generalSigma,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::sigmaModels::generalSigma::probablity() const
{
    return Foam::pow((1.0 - Foam::pow((3.764 - Zii_)/2.0, 2.5)), 0.4);
}


Foam::tmp<Foam::volScalarField> Foam::sigmaModels::generalSigma::sigma() const
{
    tmp<volScalarField> tsigma =
      volScalarField::New
      (
          "sigma0",
          mesh_,
          dimensionedScalar(dimCurrent*dimCurrent/dimEnergy, SMALL)
      );

    scalarField& sigma = tsigma.ref();

    //- temperature
    const scalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    sigma = a_ + b_*T + c_*exp(-d_/T);

    return tsigma;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaModels::generalSigma::generalSigma
(
    const fvMesh& mesh,
    const dictionary& sigmaDictionary
)
:
    sigmaModel(mesh, sigmaDictionary),
    a_(sigmaDictionary_.getOrDefault<scalar>("a", 0.0)),
    b_(sigmaDictionary_.getOrDefault<scalar>("b", 0.0)),
    c_(sigmaDictionary_.getOrDefault<scalar>("c", 0.0)),
    d_(sigmaDictionary_.getOrDefault<scalar>("d", 0.0)),
    TName_(sigmaDictionary_.getOrDefault<word>("T", "T")),
    por_(sigmaDictionary_.get<scalar>("porosity")),
    phi_(sigmaDictionary_.get<scalar>("phi")),
    Zii_(sigmaDictionary_.get<scalar>("Zii"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigmaModels::generalSigma::~generalSigma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaModels::generalSigma::correct
(
    volScalarField& sigmaField
) const
{
    const scalarField sigma = this->sigma()();
    const scalar prob = probablity();

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];

            sigmaField[cellI] = sigma[cellI]*(1.0 - por_)*phi_*prob;
        }
    }
}

// ************************************************************************* //
