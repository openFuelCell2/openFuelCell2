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
    public regIOobject
{
    //- Private Membrane Functions

    //- Disallow default bitwise copy construction
    dissolvedModel(const dissolvedModel&);

    //- Disallow default bitwise assignment
    void operator=(const dissolvedModel&);

protected:

    //- Mesh
    const fvMesh& mesh_;

    //- Lambda
    volScalarField lambda_;

    //- Water content transfer source/sink
    volScalarField dmdt_;

    //- Diffusion coefficient
    volScalarField Dwm_;

    //- Water activity
    volScalarField act_;


public:

    //- Runtime type information
    TypeName("dissolvedModel");

    //- Default name of the water name
    static const word modelName;

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
        inline const fvMesh& mesh() const;

        inline const volScalarField& lambda() const;
        inline volScalarField& lambda();

        inline const volScalarField& dmdt() const;
        inline volScalarField& dmdt();

        inline const volScalarField& act() const;
        inline volScalarField& act();

     //- IO

       //- Write
       virtual bool writeData(Ostream& os) const;

       //- Inherit read from regIOobject
       using regIOobject::read;

//       //- Read dictionary
//       virtual bool read(const dictionary& dict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dissolvedModelI.H"

#endif

// ************************************************************************* //
