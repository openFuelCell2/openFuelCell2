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
    Foam::noneDWModel.C

Description
    No dissolved water transfer is taken into account.
    This is the default selection if no dissolved water model is selected.

SourceFiles
    noneDWModel.C

\*---------------------------------------------------------------------------*/

#ifndef noneDWModel_H
#define noneDWModel_H

#include "dissolvedModel.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fuelCellSystem;

namespace dissolvedModels
{

/*---------------------------------------------------------------------------*\
                           Class dissolved Declaration
\*---------------------------------------------------------------------------*/

class noneDWModel
:
    public dissolvedModel
{
        //- Dictionary
        const dictionary& dict_;

public:

    //- Type name
        TypeName("none");

    // Constructors

        noneDWModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~noneDWModel();


    // Member Functions
        //- update
        virtual void update(const word& clName);

        //- correct
        virtual void correct();

        //- solve, if necessary
        virtual void solve();

        //- map
        virtual void mapToCell(fuelCellSystem&);

        //- map from cell
        virtual void mapFromCell(fuelCellSystem&);

    //- IO

        virtual bool read(const dictionary& dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace dissolvedModels
} // End namespace Foam

#endif

// ************************************************************************* //
