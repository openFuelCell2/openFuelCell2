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
    Foam::ActivationOverpotentialModel.C

SourceFiles
    ActivationOverpotentialModel.C

\*---------------------------------------------------------------------------*/

#ifndef ActivationOverpotentialModel_H
#define ActivationOverpotentialModel_H

#include "activationOverpotentialModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

template<class ThermoType> class pureMixture;
template<class ThermoType> class multiComponentMixture;

/*---------------------------------------------------------------------------*\
                           Class ActivationOverpotentialModel Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo>
class ActivationOverpotentialModel
:
    public activationOverpotentialModel
{
protected:

        //- Thermo
        const Thermo& thermo_;

    // Private member functions

        //- Get a reference to the local thermo for a pure mixture
        template<class ThermoType>
        const typename pureMixture<ThermoType>::thermoType&
        getLocalThermo
        (
            const word& speciesName,
            const pureMixture<ThermoType>& globalThermo
        ) const;

        //- Get a reference to the local thermo for a multi component mixture
        template<class ThermoType>
        const typename multiComponentMixture<ThermoType>::thermoType&
        getLocalThermo
        (
            const word& speciesName,
            const multiComponentMixture<ThermoType>& globalThermo
        ) const;

public:

    // Constructors

        ActivationOverpotentialModel
        (
            const phaseModel& phase,
            const dictionary& dict
        );

    //- Destructor
        virtual ~ActivationOverpotentialModel();


    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiation for multi-component-multi-component pairs
#define makeActivationOverpotentialType(Type, Thermo, Comp, Mix, Phys)\
                                                                               \
    typedef Thermo<Comp, SpecieMixture<Mix<Phys>>>                             \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    addActivationOverpotentialToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys                                         \
    )

// Addition to the run-time selection table
#define addActivationOverpotentialToRunTimeSelectionTable(Type, Thermo)\
                                                                               \
    typedef Type<Thermo>                                          \
        Type##Thermo;                                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Type##Thermo,                                             \
        (                                                                      \
            word(Type##Thermo::typeName_()) + "<"                 \
          + word(Thermo::typeName) + ">"                                  \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        activationOverpotentialModel,                                             \
        Type##Thermo,                                             \
        dictionary                                                             \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ActivationOverpotentialModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
