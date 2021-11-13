/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd
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
    Foam::wallBoilingModels:filmBoilingModels:::Bromley

Description
    Boiling film correlation.

    References:
    \verbatim
        A. Bromley, Heat transfer in stable film boiling,
        Chem. Eng. Prog. 58 (1950) 67–72.
    \endverbatim

SourceFiles
    Bromley.C

\*---------------------------------------------------------------------------*/

#ifndef Bromley_H
#define Bromley_H

#include "filmBoilingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class Bromley Declaration
\*---------------------------------------------------------------------------*/

class Bromley
:
    public filmBoilingModel
{

    // Private data:

        //- Coefficient for nucleation site density
        scalar Cn_;

        //- Wall emissivity
        scalar emissivity_;

        //- Characteristic lenght-scale
        scalar L_;

public:

    //- Runtime type information
    TypeName("Bromley");

    // Constructors

        //- Construct from a dictionary
        Bromley(const dictionary& dict);


    //- Destructor
    virtual ~Bromley();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> htcFilmBoil
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;


        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmBoilingModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
