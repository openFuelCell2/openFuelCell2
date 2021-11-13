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
    Foam::fluid.C

SourceFiles
    fluid.C

\*---------------------------------------------------------------------------*/

#ifndef fluid_H
#define fluid_H

#include "regionType.H"
#include "phaseSystem.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fvMesh.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fuelCellSystem;

namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class fluid Declaration
\*---------------------------------------------------------------------------*/

class fluid
:
    public regionType
{
protected:
        //- mesh, parent mesh
        const fvMesh& mesh_;

        //- phase type
        autoPtr<phaseSystem> phases_;


public:
        TypeName("fluid");

    // Constructors
        fluid
        (
            const fvMesh& mesh,
            const word& regionName
        );

    //- Destructor
    virtual ~fluid();

    // Member Functions
        virtual bool active() const
        {
            return true;
        }

        virtual void correct();

        virtual void setRDeltaT();

        virtual void solve();

        virtual void mapToCell(fuelCellSystem&);

        virtual void mapFromCell(fuelCellSystem&);

        phaseSystem& phases()
        {
            return phases_();
        }

        const phaseSystem& phases() const
        {
            return phases_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
