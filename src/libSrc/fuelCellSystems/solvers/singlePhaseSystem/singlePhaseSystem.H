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
    Foam::singlePhaseSystem

Description
    To solve the governing equations of single phase system, like SIMPLE(C).

SourceFiles
    singlePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseSystem_H
#define singlePhaseSystem_H

#include "dictionary.H"
// #include "phaseModel.H"
#include "volFields.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class singlePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseSystem
:
    public phaseSystem
{
    // Private data

        //- Phase model 1
        phaseModel& phase_;

public:
        TypeName("singlePhaseSystem");

    // Constructors

        //- Construct from fvMesh
        singlePhaseSystem(const fvMesh&);

    //- Destructor
    virtual ~singlePhaseSystem();

        //- Solve for the two-phase-fractions
        virtual void solve();

        //- Correct two-phase properties other than turbulence
        virtual void correct();

        // Access

            //- Single phase
            virtual bool isSinglePhase() const
            {
                return true;
            }

            //- Return relative flux magnitude
            virtual tmp<surfaceScalarField> phirMag() const;

            //- Return max magnitude flux
            virtual tmp<surfaceScalarField> phiMagMax() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
