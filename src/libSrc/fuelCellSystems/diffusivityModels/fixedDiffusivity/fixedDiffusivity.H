/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::diffusivityModels::fixedDiffusivity

Description
    A fixed value diffusivity model

    The input dictionary, diffsivity, is of the following form:

    @verbatim
    diffusivity
    {
        type        fixedDiffusivity;
        diff0       diff0 [0 2 -1 0 0 0 0]  2.0e-4;    [m^2/s]
    }
    @endverbatim


SourceFiles
    fixedDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef fixedDiffusivity_H
#define fixedDiffusivity_H

#include "diffusivityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

//class diffusivityModel;

/*---------------------------------------------------------------------------*\
                           Class fixedDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class fixedDiffusivity
:
    public diffusivityModel
{


protected:

    // protected data

        //- the fixed diffusivity
        const dimensionedScalar diff0_;


public:

    //- Runtime type information
    TypeName("fixedDiffusivity");


    // Constructors

        //- construct from mesh, diffusivity field, cell list, and dictionary
        fixedDiffusivity
        (
            const word& name,
            const fvMesh& mesh,
            scalarField& diff,
            const dictionary& dict
        );


    // Destructor

        ~fixedDiffusivity()
        {}


    // Member functions

        //- return the fixed diffusivity
        const dimensionedScalar& diff0() const
        {
            return diff0_;
        }

        //- provide feedback for user
        void writeData();

        //- evaluate the diffusivity
        void evaluate();

        //- model info
        bool isFixed()
        {
            return 1;
        }
        bool isBinary()
        {
            return 0;
        }
        bool isKnudsen()
        {
            return 0;
        }

        void setSpecies(word spA, word spB)
        {
            // not implemented
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

