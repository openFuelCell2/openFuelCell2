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
    Foam::masterRegionToCell

Description
    A topoSetSource to select zone based on masterCellZone.

SourceFiles
    masterRegionToCell.C

\*---------------------------------------------------------------------------*/

#ifndef masterRegionToCell_H
#define masterRegionToCell_H

#include "topoSetSource.H"
#include "wordRe.H"
#include "IOdictionary.H"
#include "labelIOList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class masterRegionToCell Declaration
\*---------------------------------------------------------------------------*/

class masterRegionToCell
:
    public topoSetSource
{

    // Private data

        //- Map to the master mesh
        labelIOList cellMapIO_;

        //- Map from the master mesh
        Map<label> cellMap_;

        //- Name/regular expression of cellZone
        wordRe zoneName_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("masterRegionToCell");

    // Constructors

        //- Construct from components
        masterRegionToCell
        (
            const polyMesh& mesh,
            const word& zoneName
        );

        //- Construct from dictionary
        masterRegionToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~masterRegionToCell();


    // Member Functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
