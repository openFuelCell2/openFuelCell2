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
    Foam::ChangJaffe.C

Description
    The linear relation in electrochemistry refers to a simplified
    approximation of the kinetics of electrochemical reactions.
    Unlike the more detailed Butler-Volmer and Tafel relations,
    the linear relation assumes a simplified, linear relationship
    between the rate of an electrochemical reaction and the overpotential
    (the difference between the actual electrode potential and the equilibrium potential).

    \f[
        j = j_0 \prod_i\left(\frac{p_i}{p_{ref}}\right)^\ksi (1 - s)^\gamma \eta
    \f]

    Read more:
    \vartable
        Chang, Hung‐Chi, and George Jaffé. "Polarization in electrolytic solutions.
        Part I. Theory." The Journal of Chemical Physics 20.7 (1952): 1071-1077.
    \endvartable

SourceFiles
    ChangJaffe.C

\*---------------------------------------------------------------------------*/

#ifndef ChangJaffe_H
#define ChangJaffe_H

#include "ActivationOverpotentialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

namespace activationOverpotentialModels
{

/*---------------------------------------------------------------------------*\
                           Class ChangJaffe Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo>
class ChangJaffe
:
    public ActivationOverpotentialModel<Thermo>
{

public:
    //- Runtime type information
    TypeName("ChangJaffe");

    // Constructors

    ChangJaffe
    (
        const phaseModel& phase,
        const dictionary& dict
    );

    //- Destructor
        virtual ~ChangJaffe();


    // Member Functions

        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace ActivationOverpotentialModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ChangJaffe.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
