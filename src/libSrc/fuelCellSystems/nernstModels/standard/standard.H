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

Class
    Foam::standard.C

SourceFiles
    standard.C

\*---------------------------------------------------------------------------*/

#ifndef standard_H
#define standard_H

#include "NernstModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

namespace nernstModels
{

/*---------------------------------------------------------------------------*\
                           Class standard Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo, class OtherThermo>
class standard
:
    public NernstModel<Thermo, OtherThermo>
{
    //- Phase 1
    const phaseModel& phase1_;

    //- Phase 2
    const phaseModel& phase2_;

public:
    //- Runtime type information
    TypeName("standard");

    // Constructors

    standard
    (
        const phaseModel& phase1,
        const phaseModel& phase2,
        const dictionary& dict
    );

    //- Destructor
        virtual ~standard();


    // Member Functions

        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace NernstModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "standard.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //