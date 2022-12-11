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
    Foam::driftFluxSystem

Description

SourceFiles
    driftFluxSystem.C

\*---------------------------------------------------------------------------*/

#ifndef driftFluxSystem_H
#define driftFluxSystem_H

#include "dictionary.H"
#include "twoPhaseInteractingMixture.H"
#include "relativeVelocityModel.H"
#include "CompressibleTurbulenceModel.H"
#include "volFields.H"
#include "phaseSystem.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class driftFluxSystem Declaration
\*---------------------------------------------------------------------------*/

class driftFluxSystem
:
    public phaseSystem
{
    // Private data

        //- Phase models
        phaseModel& phase1_;
        phaseModel& phase2_;

        // Mixture model
        autoPtr<twoPhaseInteractingMixture> mixture_;

        //- Mass averaged velocity
        volVectorField& U_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture phi
        surfaceScalarField rhoPhi_;

        //- Capillary pressure
        volTensorField pc_;

        //- dpc/ds derivation of pc to liquid water saturation
        volTensorField dpcds_;

        //- RelativeVelocityModel
        autoPtr<relativeVelocityModel> UdmModel_;

        //- Turbulence model
        autoPtr
        <
            CompressibleTurbulenceModel<twoPhaseInteractingMixture>
        > turbulence_;

public:
        TypeName("driftFluxSystem");

    // Constructors

        //- Construct from fvMesh
        driftFluxSystem(const fvMesh&);

    //- Destructor
    virtual ~driftFluxSystem();

        // Member Functions

            //- Return phase model 1
            const phaseModel& phase1() const;

            //- Access phase model 1
            phaseModel& phase1();

            //- Return phase model 2
            const phaseModel& phase2() const;

            //- Access phase model 2
            phaseModel& phase2();

            //- Return the phase not given as an argument
            const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Solve for the phase fractions
            virtual void solve();

            //- is Single phase
            virtual bool isSinglePhase() const;


        // Transfers

            //- Return the heat transfer matrices
            virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

            //- Return the mass transfer matrices
            virtual autoPtr<massTransferTable> massTransfer() const = 0;

            //- Return max flux magnitude
            virtual tmp<surfaceScalarField> phiMagMax() const;

            //- Return relative magnitude flux
            virtual tmp<surfaceScalarField> phirMag() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "driftFluxSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
