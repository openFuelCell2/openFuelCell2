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

#include "lambdaSigma.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{
    defineTypeNameAndDebug(lambdaSigma, 0);

    addToRunTimeSelectionTable
    (
        sigmaModel,
        lambdaSigma,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaModels::lambdaSigma::lambdaSigma
(
    const fvMesh& mesh,
    const dictionary& sigmaDictionary
)
:
    sigmaModel(mesh, sigmaDictionary),
    lambdaName_(sigmaDictionary_.lookup("lambdaName")),
    TName_(sigmaDictionary_.lookupOrDefault<word>("T", "T")),
    catalyst_(sigmaDictionary_.lookupOrDefault<Switch>("catalyst", false)),
    corr_(sigmaDictionary_.lookupOrDefault<scalar>("corr", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigmaModels::lambdaSigma::~lambdaSigma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaModels::lambdaSigma::correct
(
    volScalarField& sigmaField
) const
{
//     word lambdaName = IOobject::groupName(lambdaName_, mesh_.name());
    const scalarField& lambda = mesh_.lookupObject<volScalarField>(lambdaName_);

//     word TName = IOobject::groupName(TName_, mesh_.name());
    const scalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const scalar w = lambda[cellI];

            if(w >= 1.0)
            {
                sigmaField[cellI] = (0.514*w - 0.326)
                            *exp(1268.*(1./303. - 1./T[cellI]));
            }
            else
            {
                sigmaField[cellI] = 0.1879*w
                            *exp(1268.*(1./303. - 1./T[cellI]));
            }

            sigmaField[cellI] *= corr_;

            if(catalyst_)
            {
                scalar porNaf(sigmaDictionary_.lookupOrDefault<scalar>("porNaf", 1.0));
                scalar por(sigmaDictionary_.lookupOrDefault<scalar>("porosity", 0.0));

                sigmaField[cellI] *= pow((scalar(1) - por)*porNaf, 1.5);
            }
        }
    }
}

// ************************************************************************* //
