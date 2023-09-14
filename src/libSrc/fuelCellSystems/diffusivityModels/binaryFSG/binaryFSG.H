/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::diffusivityModels::binaryFSG

Description
    Fuller, Schettler, Giddings (1966) binary diffusion correlation

               1e-3 * T^{1.75} * sqrt(1/mA + 1/mB)
    D = 1e-4 * -----------------------------------
                  p * [ vA^{1/3} + vB^{1/3} ]^2

    where
        D = diffusivity ......... [m^2/s]
        T = temperature ......... [K]
        p = total pressure ...... [atm]
        m = molecular weight .... [kg/kmol]
        v = diffusion volume .... [cm^3]     NOTE: cm
        A,B = species index
    
    Fuller, Schettler, and Giddings,
    A new method for prediction of binary gas-phase diffusion coefficients,
    Industrial and Engineering Chemistry, v58, n5, May, 1966, pp 19-27.

    The input dictionary, diffsivity, is of the following form:

    \vartable
    diffusivity
    {
        type        binaryFSG;
        alphaName   alpha.air;        // Volume fraction, by default 1
        Tname       T;                // T (temperature) file name
        pName       p;                // p (pressure/density) file name
        speciesA    O2;               // diffusing species
        speciesB    N2;               // background species
    }
    \endvartable

    The model uses the fsgMolecularWeights and fsgDiffusionVolumes lookup
    tables to obtain the species molecular weights and diffusion volumes.

SourceFiles
    binaryFSG.C

\*---------------------------------------------------------------------------*/

#ifndef binaryFSG_H
#define binaryFSG_H

#include "diffusivityModel.H"
#include "fsgMolecularWeights.H"
#include "fsgDiffusionVolumes.H"
#include "dimensionedScalar.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

//class diffusivityModel;

/*---------------------------------------------------------------------------*\
                           Class binaryFSG Declaration
\*---------------------------------------------------------------------------*/

class binaryFSG
:
    public diffusivityModel
{


protected:

    // protected data

        //- name of temperature field
        word Tname_;

        //- name of pressure field (assume p = pressure/density)
        word pName_;

        //- species "A" name  (diffusing species)
        word spA_;

        //- species "A" molecular weight
        scalar mA_;

        //- species "A" FSG diffusion volume
        scalar vA_;

        //- species "B" name  (background species)
        word spB_;

        //- species "B" molecular weight
        scalar mB_;

        //- species "B" FSG diffusion volume
        scalar vB_;


public:

    //- Runtime type information
    TypeName("binaryFSG");


    // Constructors

        //- construct from mesh, scalarField, labelList and dictionary
        binaryFSG
        (
            const word& name,
            const fvMesh& mesh,
            scalarField& diff,
            const dictionary& dict
        );

        //- construct from components
        binaryFSG
        (
            const word& name,
            const fvMesh& mesh,
            scalarField& diff,
            word Tname_,
            word pName_,
            word spA_,
            word spB_
        );


    // Destructor

        ~binaryFSG()
        {}


    // Member functions

        //- provide feedback for user
        void writeData();

        //- evaluate the diffusivity
        void evaluate();

        //- set the species
        void setSpecies(word spA, word spB);

        //- model info
        bool isFixed()
        {
            return 0;
        }
        bool isBinary()
        {
            return 1;
        }
        bool isKnudsen()
        {
            return 0;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

