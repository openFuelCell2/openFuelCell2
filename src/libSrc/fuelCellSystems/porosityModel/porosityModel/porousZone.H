/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::porousZone

Description
    Top level model for porosity models

SourceFiles
    porousZone.C
    porousZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef porousZone_H
#define porousZone_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "dimensionedVector.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class porousZone Declaration
\*---------------------------------------------------------------------------*/

class porousZone
:
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        porousZone(const porousZone&);

        //- Disallow default bitwise assignment
        void operator=(const porousZone&);


protected:

    // Protected data

        //- Porosity name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dictionary used for model construction
        const dictionary dict_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Porosity active flag
        bool active_;

        //- Name(s) of cell-zone
        wordRe zoneName_;

        //- Cell zone IDs
        labelList cellZoneIDs_;

        //- Local co-ordinate system
        coordinateSystem coordSys_;
//Added SBB
       scalar porosity_;
//End added SBB

    // Protected Member Functions


        //- Transform the model data wrt mesh changes
        virtual void calcTranformModelData() = 0;

        //- Adjust negative resistance values to be multiplier of max value
        void adjustNegativeResistance(dimensionedVector& resist);

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const = 0;

        virtual void correct(fvVectorMatrix& UEqn) const = 0;

        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const = 0;

        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const = 0;

        //- Return label index
        label fieldIndex(const label index) const;


public:

    //- Runtime type information
    TypeName("porousZone");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        porousZone,
        mesh,
        (
            const word& modelName,
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName
        ),
        (modelName, name, mesh, dict, cellZoneName)
    );

    //- Constructor
    porousZone
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Return pointer to new porousZone object created on the freestore
    //  from an Istream
    class iNew
    {
        //- Reference to the mesh database
        const fvMesh& mesh_;
        const word& name_;

    public:

        iNew
        (
            const fvMesh& mesh,
            const word& name
        )
        :
            mesh_(mesh),
            name_(name)
        {}

        autoPtr<porousZone> operator()(Istream& is) const
        {
            const dictionary dict(is);

            return autoPtr<porousZone>
            (
                porousZone::New
                (
                    name_,
                    mesh_,
                    dict
                )
            );
        }
    };

    //- Selector
    static autoPtr<porousZone> New
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Destructor
    virtual ~porousZone();


    // Member Functions

// Added SBB
        //- Return const access to the zone name

        inline const wordRe& zoneName() const;

	//- dictionary values used for the porousZone

        inline const dictionary& dict() const;

	//- Return porosity

        inline scalar porosity() const;
// End added SBB

        //- Return const access to the porosity model name
        inline const word& name() const;

        //- Return const access to the porosity active flag
        inline bool active() const;

        //- Return const access to the cell zone IDs
//        inline const labelList& cellZoneIDs() const; Renamed SBB
        inline const labelList& zoneIds() const;//End renamed SBB

        //- Transform the model data wrt mesh changes
        virtual void transformModelData();

        //- Return the force over the cell zone(s)
        virtual tmp<vectorField> force
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        virtual void addResistance(fvVectorMatrix& UEqn);

        //- Add resistance
        virtual void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC
        );

        //- Capillary pressure . dpc/ds
        virtual void dpcds
        (
            volTensorField& pc,
            volTensorField& dpcds
        ) = 0;

        //- Update the velocities
        //- Specifically for drift-flux model
        virtual void correctU
        (
            volVectorField& U1,
            volVectorField& U2,
            const volVectorField& U,
            const volTensorField& dpcdc
        ) = 0;

    // I-O

        //- Write
        virtual bool writeData(Ostream& os) const;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Read porosity dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "porousZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
