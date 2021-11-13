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

#ifndef standardDWModel_H
#define standardDWModel_H

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

class standardDWModel
:
    public dissolvedModel
{
        //- Dictionary
        const dictionary& dict_;

        //- Water transfer coefficient
        HashTable<scalar> ksi_;

        //- nd
        dimensionedScalar nd_;

        //- rho / membrane equivillent weight (mol/m3)
        dimensionedScalar rhoOnEW_;

        //- Name of current
        word iName_;

        //- Name of temperature
        word TName_;

        //- Relax factor
        scalar relax_;

        //- Correct factor for diffusion coefficient
        scalar corr_;

        //- Functions

            //- diffusion coefficient
            scalar DNaf(const scalar&, const scalar&) const;

            //- Effective water content
            scalar effLam(const scalar&) const;


public:

    //- Runtime type information
    TypeName("standard");

    // Constructors

        standardDWModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            word modelType
        );

    //- Destructor
    virtual ~standardDWModel();


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
