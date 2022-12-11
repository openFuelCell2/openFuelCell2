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
    Foam::Constant.C

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "ActivationOverpotentialModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

namespace activationOverpotentialModels
{

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo>
class Constant
:
    public ActivationOverpotentialModel<Thermo>
{
    //- Constant values
    //- A/m3
    dimensionedScalar jj_;

public:
    //- Runtime type information
    TypeName("Constant");

    // Constructors

    Constant
    (
        const phaseModel& phase,
        const dictionary& dict
    );

    //- Destructor
        virtual ~Constant();


    // Member Functions

        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace ActivationOverpotentialModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Constant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
