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

Class
    Foam::sigmaModels::porousSigma

Description
    Electric conductivity in porous materials.
    If this is the electrode, the occupiation of Nafion should be taken into account.
    \f[
        \sigma = \sigma_0 * (1 - porosity)
    \f]
    
    otherwise:

    \f[
        \sigma = \sigma_0 * ((1 - porosity) * porNaf)^1.5
    \f]

    where
    \vartable
        \sigma_0        | Electric conductivity of electron/ion [S/m]
        porosity        | Porosity of the porous material [-]
        porNaf          | The volume fraction of Nafion material in the catalyst layer [-]
    \endvartable

SourceFiles
    porousSigma.C

\*---------------------------------------------------------------------------*/

#ifndef porousSigma_H
#define porousSigma_H

#include "sigmaModel.H"
#include "List.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{

/*---------------------------------------------------------------------------*\
                           Class porousSigma Declaration
\*---------------------------------------------------------------------------*/

class porousSigma
:
    public sigmaModel
{
    // Private data

        //- porosity
        scalar porosity_;

        //- sigma
        scalar sigma_;

        //- catalyst
        bool catalyst_;

public:

    //- Runtime type information
    TypeName("porousSigma");


    // Constructors

        //- Construct from components
        porousSigma
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    //- Destructor
    virtual ~porousSigma();


    // Member Functions

        //- Return the sigma field
        virtual void correct(volScalarField& sigmaField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sigmaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
