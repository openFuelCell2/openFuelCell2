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
    Foam::wallBoilingModels:MHFModels:::Jeschar

Description
    Minimum heat flux (MHF) model.

    References:
    \verbatim
        Jeschar, E. Specht, C. Kohler, Heat Transfer during Cooling of
        Heated Metallic Objects with Evaporating Liquids,
        Theory and Technology in Quenching, Springer, 1992. Chapter 4.
    \endverbatim

SourceFiles
    Jeschar.C

\*---------------------------------------------------------------------------*/

#ifndef Jeschar_H
#define Jeschar_H

#include "MHFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class Jeschar Declaration
\*---------------------------------------------------------------------------*/

class Jeschar
:
    public MHFModel
{

    // Private data:

        //- Burn out factor
        scalar Kmhf_;

public:

    //- Runtime type information
    TypeName("Jeschar");

    // Constructors

        //- Construct from a dictionary
        Jeschar(const dictionary& dict);


    //- Destructor
    virtual ~Jeschar();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> MHF
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
