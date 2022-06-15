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
    Foam::diffusivityModelList

Description
    List container for porosity models

SourceFiles
    diffusivityModelList.C

\*---------------------------------------------------------------------------*/

#ifndef diffusivityModelList_H
#define diffusivityModelList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class diffusivityModelList Declaration
\*---------------------------------------------------------------------------*/

class diffusivityModelList
:
    public PtrList<diffusivityModel>
{
        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- A dictionary
        const dictionary dict_;

        //- diff
        scalarField diff_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        diffusivityModelList(const diffusivityModelList&);

        //- Disallow default bitwise assignment
        void operator=(const diffusivityModelList&);


public:

    //- Constructor
    diffusivityModelList(const fvMesh& mesh, word phaseName);

    //- Constructor
    diffusivityModelList(const fvMesh& mesh, const dictionary& dict, word phaseName);

    //- Destructor
    ~diffusivityModelList();


    // Member Functions

        static const word dictName;

        //- Return active status
        bool active() const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Correct
        void correct
        (
            const PtrList<volScalarField>& X,
            PtrList<volScalarField>& diffSp,
            label inertIndex
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
