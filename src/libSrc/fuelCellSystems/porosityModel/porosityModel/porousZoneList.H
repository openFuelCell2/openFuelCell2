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
    Foam::porousZoneList

Description
    List container for porosity models

SourceFiles
    porousZoneList.C

\*---------------------------------------------------------------------------*/

#ifndef porousZoneList_H
#define porousZoneList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "porousZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class porousZoneList;
Ostream& operator<<(Ostream& os, const porousZoneList& models);

/*---------------------------------------------------------------------------*\
                      Class porousZoneList Declaration
\*---------------------------------------------------------------------------*/

class porousZoneList
:
    public PtrList<porousZone>
{

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Disallow default bitwise copy construct
        porousZoneList(const porousZoneList&);

        //- Disallow default bitwise assignment
        void operator=(const porousZoneList&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        const IOdictionary dict_;


public:

    //- Constructor
    porousZoneList(const fvMesh& mesh);

    //- Destructor
    ~porousZoneList();


    // Member Functions

        static const word dictName;

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Add resistance
        void addResistance(fvVectorMatrix& UEqn);

        //- Add resistance
        void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        );

        void dpcds
        (
            volTensorField& pc,
            volTensorField& dpcds
        );

        void correctU
        (
            volVectorField& U1,
            volVectorField& U2,
            const volVectorField& U,
            const volTensorField& dpcds
        );

        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const porousZoneList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
