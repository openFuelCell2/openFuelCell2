/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
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
    Foam::regionType.C

SourceFiles
    regionType.C

\*---------------------------------------------------------------------------*/

#ifndef regionType_H
#define regionType_H

#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Map.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fuelCellSystem;

/*---------------------------------------------------------------------------*\
                           Class regionType Declaration
\*---------------------------------------------------------------------------*/

class regionType
:
    public fvMesh
{
    // Private functions

      //- IOobject
      IOobject createIOobject(const fvMesh&) const;

protected:
        //- mesh
        const fvMesh& mesh_;

        const IOdictionary dict_;

        //- read from files
        labelIOList faceRegionAddressingIO_;

        labelIOList cellMapIO_;

        labelIOList patchesMapIO_;

        labelList faceMap_;

        scalarField faceMask_;

        //- HashTable with label and label
        Map<label> faceRegionAddressing_;

        Map<label> cellMap_;

        Map<label> patchesMap_;

public:

    //- Runtime type information
    TypeName("regionType");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            regionType,
            dictionary,
            (
                const fvMesh& mesh,
                const word& name
            ),
            (mesh, name)
        );

    // Constructors
        regionType
        (
            const fvMesh& mesh,
            const word& name
        );

    //- Destructor
    virtual ~regionType();

        // Selectors, keep this...
        static autoPtr<regionType> New
        (
            const fvMesh& mesh,
            const word& name,
            const word& modelType
        );

    // Member Functions

        //- Dictionary name
        static const word dictName;

        //- This dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        const labelList& faceRegionAddressingIO() const
        {
            return faceRegionAddressingIO_;
        }

        const Map<label>& faceRegionAddressing() const
        {
            return faceRegionAddressing_;
        }

        const labelList& faceMap() const
        {
            return faceMap_;
        }

        const scalarField& faceMask() const
        {
            return faceMask_;
        }

        const labelList& cellMapIO() const
        {
            return cellMapIO_;
        }

        const Map<label>& cellMap() const
        {
            return cellMap_;
        }

        const labelList& patchesMapIO() const
        {
            return patchesMapIO_;
        }

        const Map<label>& patchesMap() const
        {
            return patchesMap_;
        }

        // virtual functions

        virtual bool active() const
        {
            return true;
        }

        virtual void correct() = 0;

        virtual void setRDeltaT() = 0;

        virtual void solve() = 0;

        virtual void mapToCell(fuelCellSystem&) = 0;

        virtual void mapFromCell(fuelCellSystem&) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
