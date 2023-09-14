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
    Foam::electroChemicalReaction

Description
    Class for electroChemicalReactions.

TODO
    Multiple electrochemical reactions
    Combined chemical reactions

SourceFiles
    electroChemicalReaction.C

\*---------------------------------------------------------------------------*/

#ifndef electroChemicalReaction_H
#define electroChemicalReaction_H

#include "ThermoCombustion.H"
#include "saturationModel.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class activationOverpotentialModel;

namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     Class electroChemicalReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class electroChemicalReaction
:
    public ThermoCombustion<ReactionThermo>
{
    // Private Member Functions

        //- Disallow copy construct
        electroChemicalReaction(const electroChemicalReaction&);

        //- Disallow default bitwise assignment
        void operator=(const electroChemicalReaction&);

protected:

        //- Thermo model
        const ReactionThermo& thermo_;

        //- Activation overpotential model
        autoPtr<activationOverpotentialModel> eta_;

        //- Saturation model
        autoPtr<saturationModel> saturation_;

        //- Dissolved on/off
        Switch dissolved_;

public:
    //- Runtime type information
    TypeName("electroChemicalReaction");

    // Constructors

        //- Construct from components
        electroChemicalReaction
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& electroChemicalReactionProperties
        );


    //- Destructor
    virtual ~electroChemicalReaction();


    // Member Functions

        //- Correct electroChemicalReaction rate
        virtual void correct();

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "electroChemicalReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
