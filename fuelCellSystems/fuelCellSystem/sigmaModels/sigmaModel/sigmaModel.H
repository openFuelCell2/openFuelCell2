/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::sigmaModel

Description
    A2stract base-class for dispersed-phase particle diameter models.

SourceFiles
    sigmaModel.C
    newDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaModel_H
#define sigmaModel_H

#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "dimensionedTensor.H"
#include "runTimeSelectionTables.H"
#include "keyType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sigmaModel Declaration
\*---------------------------------------------------------------------------*/

class sigmaModel
{

protected:

    // Protected data
        const fvMesh& mesh_;

        const dictionary& sigmaDictionary_;

        //- Name(s) of cell-zone
        keyType zoneName_;

        //- Cell zone IDs
        labelList cellZoneIDs_;


public:

    //- Runtime type information
    TypeName("sigmaModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            sigmaModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& sigmaDictionary
            ),
            (mesh, sigmaDictionary)
        );


    // Constructors

        sigmaModel
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    //- Destructor
    virtual ~sigmaModel();


    // Selectors

        static autoPtr<sigmaModel> New
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    // Member Functions

        bool active() const
        {
            return true;
        }

        //- Return the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the dictionary
        const dictionary& sigmaDictionary() const
        {
            return sigmaDictionary_;
        }

        //- Correct the field
        virtual void correct(volScalarField&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
