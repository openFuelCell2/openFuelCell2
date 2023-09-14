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
    Foam::diffusivityModels::porousFSG

Description
    Porous medium diffusivity model consisting of
    harmonic average of effective binaryFSG and knudsen diffusivities:

                      porosity                  1
    D_{porousFSG} =  ---------- * ---------------------------
                     tortuosity   (1/DbinaryFSG + 1/Dknudsen)

    The input dictionary, diffsivity, is of the following form:

    \vartable
    diffusivity
    {
        type        porousFSG;
        Tname       T;                             // T (temperature) file name
        pName       p;                             // p (pressure/density) file
        speciesA    O2;                            // diffusing species
        speciesB    N2;                            // background species
        porosity    0.5;
        tortuosity  3.0;
        dPore       dPore [0 1 0 0 0 0 0] 8e-06;   // pore size
    }
    \endvartable

SourceFiles
    porousFSG.C

\*---------------------------------------------------------------------------*/

#ifndef porousFSG_H
#define porousFSG_H

#include "diffusivityModel.H"
#include "fsgMolecularWeights.H"
#include "fsgDiffusionVolumes.H"
#include "dimensionedScalar.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "Switch.H"
#include "binaryFSG.H"
#include "knudsen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

//class diffusivityModel;

/*---------------------------------------------------------------------------*\
                           Class porousFSG Declaration
\*---------------------------------------------------------------------------*/

class porousFSG
:
    public diffusivityModel
{


protected:

    // protected data

        //- name of temperature field
        word Tname_;

        //- name of pressure field (assume p = pressure/density)
        word pName_;

        //- name of alpha
        word alphaName_;

        //- species "A" name
        word spA_;

        //- species "A" molecular weight
        scalar mA_;

        //- species "A" FSG diffusion volume
        scalar vA_;

        //- species "B" name
        word spB_;

        //- species "B" molecular weight
        scalar mB_;

        //- species "B" FSG diffusion volume
        scalar vB_;

        //- porosity
        scalar eps_;

        //- tortuosity
        scalar tau_;

        //- pore size
        const dimensionedScalar dPore_;

        //- switch to calculate binary or not
        //  binary might be already calculated
        //Switch doBinary_;   //must be true for multiple species ...
                              //... removed 2012july23 to prevent problems


public:

    //- Runtime type information
    TypeName("porousFSG");


    // Constructors

        //- construct from mesh, diffusivity field, cell list, and dictionary
        porousFSG
        (
            const word& name,
            const fvMesh& mesh,
            scalarField& diff,
            const dictionary& dict
        );


    // Destructor

        ~porousFSG()
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
            return 1;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

