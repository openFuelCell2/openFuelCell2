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
    Foam::fuelCellSystem

Description

SourceFiles
    fuelCellSystem.C

\*---------------------------------------------------------------------------*/

#ifndef fuelCellSystem_H
#define fuelCellSystem_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionTypeList;

/*---------------------------------------------------------------------------*\
                      Class fuelCellSystem Declaration
\*---------------------------------------------------------------------------*/

class fuelCellSystem
:
    public IOdictionary
{
    // Private functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

    // Private members

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- temperature
        volScalarField T_;

        //- temperature source term
        fvScalarMatrix Qdot_;

        //- thermal conductivity
        volScalarField k_;

        //- contErr, fvm::Sp(contErr*Cp, T)
        volScalarField contErrCp_;

        //- specific heat
        volScalarField Cp_;

        //- rho*Cp
        volScalarField rhoCp_;

        //- heat flux
        surfaceScalarField rhoCpPhi_;

        //- flux
        surfaceScalarField phi_;

        //- lists of fluid, solid, electric
        autoPtr<regionTypeList> regions_;

public:

    // Constructors

        //- Construct from fvMesh
        fuelCellSystem(const fvMesh&);


    //- Destructor
    virtual ~fuelCellSystem();


    // Member Functions

        //- Correct
        void correct();

        void solve();

        void setRDeltaT();

        void mapToCell();

        void mapFromCell();

        // Access
            const volScalarField& T() const
            {
                return T_;
            }

            volScalarField& T()
            {
                return T_;
            }

            volScalarField& k()
            {
                return k_;
            }

            volScalarField& Cp()
            {
                return Cp_;
            }

            volScalarField& contErrCp()
            {
                return contErrCp_;
            }

            volScalarField& rhoCp()
            {
                return rhoCp_;
            }

            surfaceScalarField& rhoCpPhi()
            {
                return rhoCpPhi_;
            }

            surfaceScalarField& phi()
            {
                return phi_;
            }

            fvScalarMatrix& Qdot()
            {
                return Qdot_;
            }

            regionTypeList& regions();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
