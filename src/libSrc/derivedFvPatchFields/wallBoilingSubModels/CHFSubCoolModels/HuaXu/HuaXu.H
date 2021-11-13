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
    Foam::wallBoilingModels:CHFModels:::HuaXu

Description

    Critical heat flux for soob cool boiling flows.

    References:
    \verbatim
        T.C. Hua, J.J. Xu, Quenching boiling in subcooled liquid nitrogen
        for solidification of aqueous materials, Mater.
        Sci. Eng. A 292 (2000) 169–172.
    \endverbatim

SourceFiles
    HuaXu.C

\*---------------------------------------------------------------------------*/

#ifndef HuaXu_H
#define HuaXu_H

#include "CHFSubCoolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class HuaXu Declaration
\*---------------------------------------------------------------------------*/

class HuaXu
:
    public CHFSubCoolModel
{

    // Private data:

        //- Burn out factor
        scalar Kburn_;

public:

    //- Runtime type information
    TypeName("HuaXu");

    // Constructors

        //- Construct from a dictionary
        HuaXu(const dictionary& dict);


    //- Destructor
    virtual ~HuaXu();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> CHFSubCool
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
