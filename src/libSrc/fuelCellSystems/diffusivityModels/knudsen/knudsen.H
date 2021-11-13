/*---------------------------------------------------------------------------*\

Class
    Foam::diffusivityModels::knudsen

Description
    A knudsen diffusivity model

    //  D_{knudsen} = (poreDiameter/2)*97*sqrt(T/MW)  [m^2/s]
    //  where
    //      poreDiameter = [m]
    //      T ............ [K]
    //      MW ........... [kg/kmol]
    //  Geankoplis, Christie J, Transport Processes and Unit Operations,
    //  second edition (1983), Allyn and Bacon Series in Engineering,
    //  ISBN 0-205-07788-9, page 452.

    The input dictionary, diffusivity, is of the following form:

    @verbatim
    diffusivity
    {
        type     knudsen;
        Tname    T;                              // T (temperature) file name
        dPore    dPore [0 3 0 0 0 0 0] 8e-06;    // density
        MW       MW [1 0 0 0 -1 0 0] 31.9988;    // O2 [kg/kmol]
    }
    @endverbatim

SourceFiles
    knudsen.C

\*---------------------------------------------------------------------------*/

#ifndef knudsen_H
#define knudsen_H

#include "diffusivityModel.H"
#include "dimensionedScalar.H"
#include "primitiveFieldsFwd.H"
#include "labelList.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

//class diffusivityModel;

/*---------------------------------------------------------------------------*\
                           Class knudsen Declaration
\*---------------------------------------------------------------------------*/

class knudsen
:
    public diffusivityModel
{


protected:

    // protected data

        //- name of temperature field
        word Tname_;

        //- the pore diameter
        const dimensionedScalar dPore_;

        //- the molecular weight
        const dimensionedScalar MW_;


public:

    //- Runtime type information
    TypeName("knudsen");


    // Constructors

        //- construct from mesh, diffsivity field labelList and dictionary
        knudsen
        (
            const fvMesh& mesh,
            scalarField& diff,
            const labelList& cells,
            const dictionary& dict
        );

        //- construct from components
        knudsen
        (
            const fvMesh& mesh,
            scalarField& diff,
            const labelList& cells,
            word Tname,
            const dimensionedScalar& dPore,
            const dimensionedScalar& MW
        );


    // Destructor

        ~knudsen()
        {}


    // Member functions

        //- return the molecular weight
        const dimensionedScalar MW() const
        {
            return MW_;
        }

        //- return the pore diameter
        const dimensionedScalar dPore() const
        {
            return dPore_;
        }

        //- provide feedback for user
        void writeData();

        //- evaluate the diffusivity
        void evaluate();

        //- model info
        bool isFixed()
        {
            return 0;
        }
        bool isBinary()
        {
            return 0;
        }
        bool isKnudsen()
        {
            return 1;
        }

        void setSpecies(word spA, word spB)
        {
            // not implemented
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

