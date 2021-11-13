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
    Foam::wallBoilingModels:LeidenfrostModels:::Spiegler

Description
    Leidenfrost temperature model.

    References:
    \verbatim
        SPIEGLER P., HOPENFELD J., SILBERBERG M., BUMPUS J. and NORMAN A.,
        Onset of stable film boiling and the foam limit, International
        Journal of Heat and Mass Transfer, 6,11, pp.987-989, 1963
    \endverbatim

SourceFiles
    Spiegler.C

\*---------------------------------------------------------------------------*/

#ifndef Spiegler_H
#define Spiegler_H

#include "LeidenfrostModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace LeidenfrostModels
{

/*---------------------------------------------------------------------------*\
                         Class Spiegler Declaration
\*---------------------------------------------------------------------------*/

class Spiegler
:
    public LeidenfrostModel
{

private:

    //- Critical temperature
    scalar Tcrit_;


public:

    //- Runtime type information
    TypeName("Spiegler");

    // Constructors

        //- Construct from a dictionary
        Spiegler(const dictionary& dict);


    //- Destructor
    virtual ~Spiegler();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> TLeid
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

} // End namespace LeidenfrostModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
