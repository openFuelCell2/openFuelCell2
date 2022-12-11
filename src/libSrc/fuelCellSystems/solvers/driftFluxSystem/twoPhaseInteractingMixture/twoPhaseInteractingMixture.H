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
    Foam::twoPhaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    twoPhaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseInteractingMixture_H
#define twoPhaseInteractingMixture_H

#include "mixtureViscosityModel.H"
#include "compressibleTransportModel.H"
#include "phaseSystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class twoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseInteractingMixture
:
    public compressibleTransportModel
{
private:

    // Private function

        tmp<surfaceScalarField> phi(const volVectorField& U) const;

protected:

    // Protected data

        const phaseSystem& fluid_;

        word phasecName_;
        word phasedName_;

        const phaseModel& phase1_;
        const phaseModel& phase2_;

        volVectorField U_;
        surfaceScalarField phi_;

        volScalarField mu_;

public:

    TypeName("twoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        twoPhaseInteractingMixture
        (
            const phaseSystem& fluid
        );


    //- Destructor
    virtual ~twoPhaseInteractingMixture()
    {}


    // Member Functions

        //- Return continuous phase
        const phaseModel& phased() const
        {
            return fluid_.phases()[phasedName_];
        }

        phaseModel& phasedRef()
        {
            return const_cast<phaseModel&>(fluid_.phases()[phasedName_]);
        }
        const phaseModel& phasec() const
        {
            return fluid_.phases()[phasecName_];
        }

        phaseModel& phasecRef()
        {
            return const_cast<phaseModel&>(fluid_.phases()[phasecName_]);
        }

        const phaseModel& phase1() const
        {
            return phase1_;
        }

        phaseModel& phase1Ref()
        {
            return const_cast<phaseModel&>(phase1_);
        }

        const phaseModel& phase2() const
        {
            return phase2_;
        }

        phaseModel& phase2Ref()
        {
            return const_cast<phaseModel&>(phase2_);
        }

        //- Return const-access to the dispersed-phase density
        const volScalarField& rhod() const
        {
            return phased().thermo().rho();
        }

        //- Return const-access to continuous-phase density
        const volScalarField& rhoc() const
        {
            return phasec().thermo().rho();
        };

        //- Return the diameter of the dispersed-phase particles
        tmp<volScalarField> dd() const
        {
            return phased().d();
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return phased().alphaMax();
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return access to the mixture velocity
        volVectorField& URef()
        {
            return U_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const
        {
            return phasec()*rhoc() + phased()*rhod();
        }

        //- Return the mixture density for patch
        virtual tmp<scalarField> rho(const label patchi) const
        {
            return
                phased().boundaryField()[patchi]*rhod()
              + phasec().boundaryField()[patchi]*rhoc();
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            return mu_/rho();
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return mu_.boundaryField()[patchi]/rho(patchi);
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            mu_ = phasec()*phasec().thermo().mu()
                + phased()*phased().thermo().mu();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
