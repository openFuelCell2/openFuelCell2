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
    Foam::sigmaModels::generalSigma

Description
    Constant general type sigma model.
    sigma = a + b*T + c*exp(-d/T)

SourceFiles
    generalSigma.C

\*---------------------------------------------------------------------------*/

#ifndef generalSigma_H
#define generalSigma_H

#include "sigmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{

/*---------------------------------------------------------------------------*\
                           Class generalSigma Declaration
\*---------------------------------------------------------------------------*/

class generalSigma
:
    public sigmaModel
{
    // Private data

        //- a, b, c, d
        //- sigma = a + b*T + c*exp(-d/T)
        scalar a_;
        scalar b_;
        scalar c_;
        scalar d_;

        //- temperature name
        word TName_;

        //- Porosity
        scalar por_;

        //- Volume fraction
        scalar phi_;

        //- Coordinate number
        scalar Zii_;

    // Private function

        //- P value
        scalar probablity() const;

        //- Conductivity (S/m)
        tmp<volScalarField> sigma() const;


public:

    //- Runtime type information
    TypeName("generalSigma");


    // Constructors

        //- Construct from components
        generalSigma
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    //- Destructor
    virtual ~generalSigma();


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
