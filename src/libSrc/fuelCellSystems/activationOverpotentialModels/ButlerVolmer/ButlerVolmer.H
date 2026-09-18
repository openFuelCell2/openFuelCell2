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
    Foam::ButlerVolmer.C

Description
    The Butler-Volmer relation is an abstract concept in electrochemistry
    that describes the mathematical relationship between the rate of
    an electrochemical reaction and the overpotential
    (the deviation from the equilibrium potential) experienced
    by an electrode during that reaction.

    \f[
        j = j_0 \prod_i\left(\frac{p_i}{p_{ref}}\right)^\ksi (1-s)^\gamma
        \left(
            \exp{\frac{n\alpha F}{RT}\eta}
          - \exp{-\frac{n(1 - \alpha)F}{RT}\eta}
        \right)
    \f]

SourceFiles
    ButlerVolmer.C

\*---------------------------------------------------------------------------*/

#ifndef ButlerVolmer_H
#define ButlerVolmer_H

#include "ActivationOverpotentialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

namespace activationOverpotentialModels
{

/*---------------------------------------------------------------------------*\
                           Class ButlerVolmer Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo>
class ButlerVolmer
:
    public ActivationOverpotentialModel<Thermo>
{
    // Private data

        //- Upper bound on the magnitude of the exponent argument
        //  n alpha F eta / (R T) applied to both the anodic and cathodic
        //  Boltzmann factors. exp(709) is the double-precision overflow
        //  limit, so an unbounded argument turns a transiently large
        //  overpotential into inf and then into a hardware floating point
        //  exception. The bound is a numerical safety net: it is inactive
        //  whenever |eta| stays below ~maxExpArg R T / (n alpha F).
        scalar maxExpArg_;

        //- Upper bound on the magnitude of the activation overpotential (V).
        //  Bounding the exponent alone is not sufficient: the activation heat
        //  source posted to the energy equation is j*eta, and an unconverged
        //  potential field can make eta enormous while the clamped j stays
        //  finite, so j*eta still overflows and destroys the temperature field.
        //  A real cell never exceeds a few volts of activation overpotential, so
        //  this bound is inactive for any converged solution.
        scalar maxEta_;


public:
    //- Runtime type information
    TypeName("ButlerVolmer");

    // Constructors

    ButlerVolmer
    (
        const phaseModel& phase,
        const dictionary& dict
    );

    //- Destructor
        virtual ~ButlerVolmer();


    // Member Functions

        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace ActivationOverpotentialModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ButlerVolmer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
