/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::sigmaModelList

Description
    List container for sigma models

SourceFiles
    sigmaModelList.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaModelList_H
#define sigmaModelList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "sigmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class sigmaModelList;

/*---------------------------------------------------------------------------*\
                      Class sigmaModelList Declaration
\*---------------------------------------------------------------------------*/

class sigmaModelList
:
    public PtrList<sigmaModel>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sigmaModelList(const sigmaModelList&);

        //- Disallow default bitwise assignment
        void operator=(const sigmaModelList&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        const dictionary dict_;


public:

    //- Constructor
    sigmaModelList(const fvMesh& mesh);

    //- Constructor
    sigmaModelList(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    ~sigmaModelList();

    // Member Functions

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        void correct(volScalarField& sigmaField);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
