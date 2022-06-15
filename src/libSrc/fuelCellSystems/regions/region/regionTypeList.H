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
    Foam::regionTypeList

Description
    List container for regions

SourceFiles
    regionTypeList.C

\*---------------------------------------------------------------------------*/

#ifndef regionTypeList_H
#define regionTypeList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "regionProperties.H"
#include "regionType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class fuelCellSystem;

/*---------------------------------------------------------------------------*\
                      Class regionTypeList Declaration
\*---------------------------------------------------------------------------*/

class regionTypeList
:
    public PtrList<regionType>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionTypeList(const regionTypeList&);

        //- Disallow default bitwise assignment
        void operator=(const regionTypeList&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        regionProperties region_;


public:

    //- Constructor
    regionTypeList(const fvMesh& mesh);

    //- Destructor
    ~regionTypeList();

    // Member Functions

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const regionProperties& region);

        void correct();

        void setRDeltaT();

        void solve();

        void mapToCell(fuelCellSystem&);

        void mapFromCell(fuelCellSystem&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
