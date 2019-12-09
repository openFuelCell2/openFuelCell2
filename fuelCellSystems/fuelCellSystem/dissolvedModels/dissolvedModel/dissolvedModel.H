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
    Foam::epRegionModel.C

SourceFiles
    epRegionModel.C

\*---------------------------------------------------------------------------*/

#ifndef dissolvedModel_H
#define dissolvedModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fuelCellSystem;

/*---------------------------------------------------------------------------*\
                           Class dissolvedModel Declaration
\*---------------------------------------------------------------------------*/

class dissolvedModel
:
    public volScalarField
{
protected:

    //- Mesh
    const fvMesh& mesh_;

    //- Water content transfer source/sink
    volScalarField dmdt_;

    //- Diffusion coefficient
    volScalarField Dwm_;

    //- Water activity
    volScalarField act_;


public:

    //- Runtime type information
    TypeName("dissolvedModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dissolvedModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                word modelType
            ),
            (mesh, dict, modelType)
        );

    // Constructors
        dissolvedModel
        (
            const fvMesh& mesh
        );

        // Selectors, keep this...
        static autoPtr<dissolvedModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~dissolvedModel();


    // Member Functions
        //- update
        virtual void update(const word& clName) = 0;

        //- correct
        virtual void correct() = 0;

        //- solve, if necessary
        virtual void solve() = 0;

        //- map
        virtual void mapToCell(fuelCellSystem&) = 0;

        //- map from cell
        virtual void mapFromCell(fuelCellSystem&) = 0;

    // Access

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        const volScalarField& dmdt() const
        {
            return dmdt_; 
        }

        volScalarField& dmdt()
        {
            return dmdt_;
        }

        const volScalarField& act() const
        {
            return act_;
        }

        volScalarField& act()
        {
            return act_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
