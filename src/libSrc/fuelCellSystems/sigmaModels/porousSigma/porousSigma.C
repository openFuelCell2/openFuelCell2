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

#include "porousSigma.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{
    defineTypeNameAndDebug(porousSigma, 0);

    addToRunTimeSelectionTable
    (
        sigmaModel,
        porousSigma,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaModels::porousSigma::porousSigma
(
    const fvMesh& mesh,
    const dictionary& sigmaDictionary
)
:
    sigmaModel(mesh, sigmaDictionary),
    catalyst_(sigmaDictionary_.lookupOrDefault<bool>("catalyst", false))
{
    sigmaDictionary_.lookup("porosity") >> porosity_;
    sigmaDictionary_.lookup("sigma") >> sigma_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigmaModels::porousSigma::~porousSigma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaModels::porousSigma::correct
(
    volScalarField& sigmaField
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];

            sigmaField[cellI] = sigma_*(1.0 - porosity_);

            if(catalyst_)
            {
                scalar porNaf(sigmaDictionary_.lookupOrDefault<scalar>("porNaf", 1.0));

                sigmaField[cellI] = sigma_
                    *pow((scalar(1) - porosity_)*porNaf, 1.5);
            }
        }
    }
}

// ************************************************************************* //
