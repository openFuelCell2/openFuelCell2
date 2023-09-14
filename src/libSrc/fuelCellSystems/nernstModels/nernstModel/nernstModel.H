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
    Foam::nernstModel.C

Description
    The nernst potential calculated on each electrode.

SourceFiles
    nernstModel.C
    newNernstModel.C

\*---------------------------------------------------------------------------*/

#ifndef nernstModel_H
#define nernstModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "HashTable.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

/*---------------------------------------------------------------------------*\
                           Class nernstModel Declaration
\*---------------------------------------------------------------------------*/
class nernstModel
:
    public volScalarField
{
protected:

        //- Dictionary for ep parameters
        const dictionary& dict_;

        //- HashTable, species name and electron transfer
        HashTable<scalar> rxnList_;

        //- delta H
        scalarField deltaH_;

        //- delta S
        scalarField deltaS_;

        //- residual Y
        scalar residualY_;

        //- Reference pressure
        dimensionedScalar pRef_;

public:
    //- Runtime type information
    ClassName("nernstModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            nernstModel,
            dictionary,
            (
                const phaseModel& phase1,
                const phaseModel& phase2,
                const dictionary& dict
            ),
            (phase1, phase2, dict)
        );

    // Constructors

        nernstModel
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const dictionary& dict
        );

    //- Destructor
        virtual ~nernstModel();

    // Selectors

        static autoPtr<nernstModel> New
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const dictionary& dict
        );

    // Member Functions

        virtual void correct() = 0;

        const dictionary& dictnernstModel() const
        {
            return dict_;
        }

        //- Reference pressure
        const dimensionedScalar& pRef() const
        {
            return pRef_;
        }

        //- J/mol
        const scalarField& deltaH() const
        {
            return deltaH_;
        }

        //- J/(mol K)
        const scalarField& deltaS() const
        {
            return deltaS_;
        }

        const HashTable<scalar> rxnList() const
        {
            return rxnList_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
