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

Description
    The 'electronIon' region considers electron and ion/proton transfer in electrochemical devices.
    Fuel cell and electrolyzer can operate under galvanostatic or potentostatic mode, which means
    the current density or voltage is prescribed, respectively.
    Possion equations are solved accordingly in electrically conductive regions.

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

        //- Current density (A/m2)
        volVectorField i_;

        //- Potential (V)
        volScalarField phi_;

        //- Source term (A/m3)
        volScalarField j_;

        //- Electric conductivity
        volScalarField sigmaField_;

        //- Temperature
        volScalarField T_;

        //- Sigma model (electric conductivity)
        autoPtr<sigmaModelList> sigma_;

        //- Relaxation
        scalar relax_;

        //- Active side?
        //- Reference side: voltage = 0.0
        //- Active side: voltage: fixed or calculated
        Switch active_;

        //- Dissolved water
        //- In PEM fuel cell/electrolyzer, water can dissolve from one side
        //- to the other. This process is usually affected by the proton transfer.
        //- Refer to the region of CCM (carbon coated membrane)
        Switch dissolveOnOff_;
        autoPtr<dissolvedModel> dissolved_;

        //- Patch name
        //- The reference boundary for the prescribed current density or voltage
        word patchName_;

        //- Operation mode: (currentCollector_ true)
        Switch galvanostatic_;

        //- Mean current density
        //- Make sense when galvanostatic_ is true
        scalar ibar_;

        //- Voltage
        //- Make semse when galvanostatic_ is false
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

        //- Correct fields
        virtual void correct();

        //- Set reversed delta t
        virtual void setRDeltaT();

        //- Solve PDEs
        virtual void solve();

        //- Map heat source/sink to the parent mesh
        virtual void mapToCell(fuelCellSystem&);

        //- Map temperature from the parent mesh
        virtual void mapFromCell(fuelCellSystem&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
