/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd
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
    Foam::wallBoilingModels::LeidenfrostModel

Description
    Base class for nucleation site density models

SourceFiles
    LeidenfrostModel.C
    newLeidenfrostModel.C

\*---------------------------------------------------------------------------*/

#ifndef LeidenfrostModel_H
#define LeidenfrostModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class LeidenfrostModel Declaration
\*---------------------------------------------------------------------------*/

class LeidenfrostModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        LeidenfrostModel(const LeidenfrostModel&);

        //- Disallow default bitwise assignment
        void operator=(const LeidenfrostModel&);


public:

    //- Runtime type information
    TypeName("LeidenfrostModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        LeidenfrostModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        LeidenfrostModel();


    // Selectors

        //- Select null constructed
        static autoPtr<LeidenfrostModel> New(const dictionary& dict);


    //- Destructor
    virtual ~LeidenfrostModel();


    // Member Functions

        //- Calculate temperature
        virtual tmp<scalarField> TLeid
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const = 0;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
