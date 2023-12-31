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
    Foam::PurePhaseModel

Description
    Class which represents pure phases, i.e. without any species. Returns an
    empty list of mass fractions.

SourceFiles
    PurePhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef PurePhaseModel_H
#define PurePhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PurePhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class PurePhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Empty mass fraction field list
        PtrList<volScalarField> Y_;

        //- Empty mass fraction field list
        PtrList<volScalarField> X_;

        //- iDmdt, mass source/sink term
        //- Only created here
        //- Update elsewhere
        PtrList<volScalarField> iDmdt_;


public:

    // Constructors

        PurePhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~PurePhaseModel();


    // Member Functions

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
    #include "PurePhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
