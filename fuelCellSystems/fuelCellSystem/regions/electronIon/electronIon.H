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
    Foam::electronIon.C

SourceFiles
    electronIon.C

\*---------------------------------------------------------------------------*/

#ifndef electronIon_H
#define electronIon_H

#include "dictionary.H"
#include "regionType.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "wordList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fuelCellSystem;
class sigmaModelList;
class dissolvedModel;

namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class electronIon Declaration
\*---------------------------------------------------------------------------*/

class electronIon
:
    public regionType
{
private:
        //-mesh
        const fvMesh& mesh_;

        //- Dictionary for ep parameters
        const dictionary& dict_;

        //- current density
        volVectorField i_;

        //- potential
        volScalarField phi_;

        //- source term
        volScalarField source_;

        //- electronIonnic conductivity
        volScalarField sigmaField_;

        //- temperature
        volScalarField T_;

        //- sigma model
        autoPtr<sigmaModelList> sigma_;

        //- coefficient
        scalar relax_;

        //- Active side?
        //- Reference side: voltage = 0.0
        //- Active side: voltage: fixed or calculated
        bool active_;

        //- dissolved water
        autoPtr<dissolvedModel> dissolved_;

        //- patch name (active_: true)
        word patchName_;

        //- Operation mode: (currentCollector_ true)
        Switch galvanostatic_;

        // mean current density
        // only make sense when galvanostatic_ is true
        scalar ibar_;

        // voltage
        // only make semse when galvanostatic_ is false
        scalar voltage_;

public:
        TypeName("electronIon");

    // Constructors

        electronIon
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
        virtual ~electronIon();


    // Member Functions
        virtual bool active() const
        {
            return active_;
        }

        virtual void correct();

        virtual void setRDeltaT();

        virtual void solve();

        virtual void mapToCell(fuelCellSystem&);

        virtual void mapFromCell(fuelCellSystem&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
