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
    Foam::wallBoilingModels:CHFModels:::Zuber

Description
    Critical heat flux (CHF) correlation

    References:
    \verbatim
        N. Zuber, On the stability of boiling heat transfer,
        Trans. ASME 80 (1958) 711
    \endverbatim

SourceFiles
    Zuber.C

\*---------------------------------------------------------------------------*/

#ifndef Zuber_H
#define Zuber_H

#include "CHFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class Zuber Declaration
\*---------------------------------------------------------------------------*/

class Zuber
:
    public CHFModel
{

    // Private data:

        //- Coefficient constant
        scalar Cn_;

public:

    //- Runtime type information
    TypeName("Zuber");

    // Constructors

        //- Construct from a dictionary
        Zuber(const dictionary& dict);


    //- Destructor
    virtual ~Zuber();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> CHF
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

} // End namespace CHFModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
