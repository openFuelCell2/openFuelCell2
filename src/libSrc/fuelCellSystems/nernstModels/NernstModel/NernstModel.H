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
    Foam::NernstModel.C

SourceFiles
    NernstModel.C

\*---------------------------------------------------------------------------*/

#ifndef NernstModel_H
#define NernstModel_H

#include "nernstModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

template<class ThermoType> class pureMixture;
template<class ThermoType> class multiComponentMixture;

/*---------------------------------------------------------------------------*\
                           Class NernstModel Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo, class OtherThermo>
class NernstModel
:
    public nernstModel
{
protected:

        //- Thermo
        const Thermo& thermo_;

        //- Other Thermo
        const OtherThermo& otherThermo_;

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

        NernstModel
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const dictionary& dict
        );

    //- Destructor
        virtual ~NernstModel();


    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiation for multi-component-multi-component pairs
#define makeNernstType(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, SpecieMixture<Mix<Phys>>>                             \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, OtherMix<OtherPhys>>                        \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addNernstToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )

// Instantiation for multi-component-single-component pairs
#define makeSpecieNernstType(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, SpecieMixture<Mix<Phys>>>                             \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, SpecieMixture<OtherMix<OtherPhys>>>         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addNernstToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )

// Addition to the run-time selection table
#define addNernstToRunTimeSelectionTable(Type, Thermo, OtherThermo)\
                                                                               \
    typedef Type<Thermo, OtherThermo>                                          \
        Type##Thermo##OtherThermo;                                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Type##Thermo##OtherThermo,                                             \
        (                                                                      \
            word(Type##Thermo##OtherThermo::typeName_()) + "<"                 \
          + word(Thermo::typeName) + ","                                       \
          + word(OtherThermo::typeName) + ">"                                  \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        nernstModel,                                             \
        Type##Thermo##OtherThermo,                                             \
        dictionary                                                             \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NernstModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
