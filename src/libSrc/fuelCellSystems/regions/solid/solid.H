/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | openFuelCell
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
    Foam::solid.C

SourceFiles
    solid.C

\*---------------------------------------------------------------------------*/

#ifndef solid_H
#define solid_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "fuelCellSystem.H"
#include "solidThermo.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class solid Declaration
\*---------------------------------------------------------------------------*/

class solid
:
    public regionType
{
        //- mesh
        const fvMesh& mesh_;

       //- solid thermo model
       autoPtr<solidThermo> thermo_;

        //- Heat source/sink
        volScalarField Qdot_;

public:
        TypeName("solid");

    // Constructors

        solid
        (
            const fvMesh& mesh,
            const word& regionName
        );

    //- Destructor
    virtual ~solid();


    // Member Functions
        // Virtual functions, for access
        virtual bool active() const
        {
            return true;
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
