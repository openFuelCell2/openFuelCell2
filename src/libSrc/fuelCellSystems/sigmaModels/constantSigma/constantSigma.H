/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::sigmaModels::constantSigma

Description
    Constant dispersed-phase particle sigma model.

SourceFiles
    constantSigma.C

\*---------------------------------------------------------------------------*/

#ifndef constantSigma_H
#define constantSigma_H

#include "sigmaModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{

/*---------------------------------------------------------------------------*\
                           Class constantSigma Declaration
\*---------------------------------------------------------------------------*/

class constantSigma
:
    public sigmaModel
{
    // Private data

        //- The constant sigma of the region
        dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("constantSigma");


    // Constructors

        //- Construct from components
        constantSigma
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    //- Destructor
    virtual ~constantSigma();


    // Member Functions

        //- Return the sigma as a field
        virtual void correct(volScalarField& sigmaField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sigmaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
