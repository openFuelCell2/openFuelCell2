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

\*---------------------------------------------------------------------------*/
#include "singlePhaseSystem.H"
#include "addToRunTimeSelectionTable.H"

#include "fvCFD.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseSystem, 0);
    addToRunTimeSelectionTable(phaseSystem, singlePhaseSystem, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseSystem::singlePhaseSystem
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),

    phase_(phaseModels_[0])
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseSystem::~singlePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::singlePhaseSystem::phirMag() const
{
    return mag(phase_.phi());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::singlePhaseSystem::phiMagMax() const
{
    return mag(phase_.phi());
}


void Foam::singlePhaseSystem::solve()
{
    Info << "\nSolve for single phase flow:" << endl;

    fvMesh& mesh = const_cast<fvMesh&>(mesh_);

    const Time& runTime = mesh.time();

    #include "./createFields.H"
    #include "./createRhoUfIfPresent.H"

    Switch Y
    (
        pimple_.dict().lookupOrDefault<Switch>("Y", true)
    );

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
//        if (pimple_.firstIter() && !pimple_.SIMPLErho())
//        {
//            #include "rhoEqn.H"
//        }

        correct();

        #include "./YEqns.H"
        #include "./UEqn.H"

        while (pimple_.correct())
        {
            if (pimple_.consistent())
            {
                #include "./pcEqn.H"
            }
            else
            {
                #include "./pEqn.H"
            }
        }

        correctKinematics();

        if (pimple_.turbCorr())
        {
            correctTurbulence();
        }
    }
}


void Foam::singlePhaseSystem::correct()
{
    phase_.correct();
}

// ************************************************************************* //
