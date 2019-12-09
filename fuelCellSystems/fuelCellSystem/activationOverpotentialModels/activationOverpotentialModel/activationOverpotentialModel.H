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
    Foam::activationOverpotentialModel.C

SourceFiles
    activationOverpotentialModel.C

\*---------------------------------------------------------------------------*/

#ifndef activationOverpotentialModel_H
#define activationOverpotentialModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "HashTable.H"
#include "volFields.H"
#include "regionType.H"
#include "phaseModel.H"
#include "nernstModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class activationOverpotentialModel Declaration
\*---------------------------------------------------------------------------*/
class activationOverpotentialModel
:
    public volScalarField
{
protected:

        //- Phase
        const phaseModel& phase_;

        //- Dictionary
        const dictionary& dict_;

        //- Phase change
        Switch phaseChange_;

        //- Nernst model
        autoPtr<nernstModel> nernst_;

        //- Acitve zone name
        word zoneName_;

        //- Catalyst layer adjacent regions
        const dictionary& regions_;

        // Species dictionary
        // e.g. O2
        // e.g. H2
        const dictionary& species_;

        // phi names for electric field
        // electron: phi
        // ion: phi
        HashTable<word> phiNames_;

        // Properties
        // relax, alpha, j0, Xref, gamma
        scalar relax_;
        scalar alpha_;
        scalar gamma_;
        dimensionedScalar j0_;

        //- Current density
        //- A/cm3
        volScalarField j_;

        //- Create nernst model
        void createNernst();

public:
    //- Runtime type information
    TypeName("activationOverpotentialModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            activationOverpotentialModel,
            dictionary,
            (
                const phaseModel& phase,
                const dictionary& dict
            ),
            (phase, dict)
        );

    // Constructors

        activationOverpotentialModel
        (
            const phaseModel& phase,
            const dictionary& dict
        );

    //- Destructor
        virtual ~activationOverpotentialModel();

    // Selectors

        static autoPtr<activationOverpotentialModel> New
        (
            const phaseModel& phase,
            const dictionary& dict
        );

    // Member Functions

        //- Lookup regionType
        const regionType& region(const word& name) const
        {
            return phase_.mesh().time().
                lookupObject<regionType>(name);
        }

        //- Dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Phase change on/off
        const Switch& phaseChange() const
        {
            return phaseChange_;
        }

        //- Nernst model
        const nernstModel& nernst() const
        {
            return nernst_();
        }

        //- Relax factor
        const scalar& relax() const
        {
            return relax_;
        }

        //- Transfer coefficient
        const scalar& alpha() const
        {
            return alpha_;
        }

        //- Regions
        const dictionary& regions() const
        {
            return regions_;
        }

        //- Name of active regions
        const word& zoneName() const
        {
            return zoneName_;
        }

        //- Exchange current density
        //- A/cm3
        const dimensionedScalar& j0() const
        {
            return j0_;
        }

        //- Current density
        //- A/cm3
        const volScalarField& j() const
        {
            return j_;
        }

        volScalarField& j()
        {
            return j_;
        }

        //- Heat source
        virtual tmp<volScalarField> Qdot() const;

        //- Update
        virtual void correct() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
