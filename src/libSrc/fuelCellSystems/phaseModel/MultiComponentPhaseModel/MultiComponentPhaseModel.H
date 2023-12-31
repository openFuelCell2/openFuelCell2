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
    Foam::MultiComponentPhaseModel

Description
    Class which represents a phase with multiple species. Returns the species'
    mass fractions, and their governing equations.

    The diffusion coefficient is described by Fick's law.

SourceFiles
    MultiComponentPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MultiComponentPhaseModel_H
#define MultiComponentPhaseModel_H

#include "phaseModel.H"
#include "diffusivityModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MultiComponentPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MultiComponentPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- diffusion coefficient
        PtrList<volScalarField> diffSp_;

        //- diffusion models
        autoPtr<diffusivityModelList> diffModels_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- Inert species index
        label inertIndex_;

        //- Pointer list to active species
        UPtrList<volScalarField> YActive_;

        //- mole fraction
        PtrList<volScalarField> X_;

        //- Mole fraction list
        //- Note: iDmdt is created here
        //- iDmdt is updated in electroChemicalReaction
        PtrList<volScalarField> iDmdt_;

        //- Pointer list to active species
        UPtrList<volScalarField> XActive_;


public:

    // Constructors

        MultiComponentPhaseModel
        (
            const phaseSystem& phase,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~MultiComponentPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi);

            //- Return the source/sink term
            virtual const PtrList<volScalarField>& iDmdt() const;

            //- Return the source/sink term by name
            virtual const volScalarField& iDmdt(const word& name) const;

            //- Access the source/sink term
            virtual PtrList<volScalarField>& iDmdtRef();

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef();

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef();

            //- Return the species mole fractions
            virtual const PtrList<volScalarField>& X() const;

            //- Return a species mole fraction by name
            virtual const volScalarField& X(const word& name) const;

            //- Access the species mole fractions
            virtual PtrList<volScalarField>& XRef();

            //- Return the active species mole fractions
            virtual const UPtrList<volScalarField>& XActive() const;

            //- Access the active species mole fractions
            virtual UPtrList<volScalarField>& XActiveRef();

            //- Return the total mass transfer
            virtual tmp<volScalarField> dmdt() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MultiComponentPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
