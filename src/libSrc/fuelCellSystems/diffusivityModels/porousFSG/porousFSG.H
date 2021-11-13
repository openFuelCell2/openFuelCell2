/*---------------------------------------------------------------------------*\

Class
    Foam::diffusivityModels::porousFSG

Description
    Porous medium diffusivity model consisting of
    harmonic average of effective binaryFSG and knudsen diffusivities:

                      porosity                  1
    D_{porousFSG} =  ---------- * ---------------------------
                     tortuosity   (1/DbinaryFSG + 1/Dknudsen)

    The input dictionary, diffsivity, is of the following form:

    @verbatim
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
    @endverbatim

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
            const fvMesh& mesh,
            scalarField& diff,
            const labelList& cells,
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

