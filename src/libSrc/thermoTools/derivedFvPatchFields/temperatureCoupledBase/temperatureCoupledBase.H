/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019,2021 OpenCFD Ltd.
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
    Foam::temperatureCoupledBase

Description
    Common functions used in temperature coupled boundaries.

    The thermal conductivity \c kappa may be obtained by the following methods:
      - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
        defined by 'kappa'
      - 'fluidThermo' : use fluidThermo and default
        compressible::turbulenceModel to calculate kappa
      - 'solidThermo' : use solidThermo kappa()
      - 'directionalSolidThermo': uses look up for volSymmTensorField for
        transformed kappa vector. Field name definable in 'alphaAni',
        named 'Anialpha' in solid solver by default
      - 'function' : kappa, alpha directly specified as Function1

    \par Keywords provided by this class:
      \table
        Property     | Description                          | Required | Default
        kappaMethod  | Thermal conductivity method          | yes |
        kappa        | Name of thermal conductivity field   | partly |
        alpha        | Name of thermal diffusivity field    | partly |
        alphaAni     | Name of non-isotropic alpha          | partly |
        kappaValue   | Function1 supplying kappa            | partly |
        alphaValue   | Function1 supplying alpha            | partly |
      \endtable

Usage
    \verbatim
    nonIsotropicWall
    {
        ...
        kappaMethod     directionalSolidThermo;
        alphaAni        Anialpha;
        ...
    }
    \endverbatim

    \verbatim
    specifiedWall
    {
        ...
        kappaMethod     function;
        kappaFunction   constant 1.0;
        alphaFunction   constant 100.0;
        ...
    }
    \endverbatim

SourceFiles
    temperatureCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureCoupledBase_H
#define temperatureCoupledBase_H

#include "scalarField.H"
#include "Enum.H"
#include "fvPatch.H"
#include "PatchFunction1.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class temperatureCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class temperatureCoupledBase
{
public:

    // Public Enumerations

        //- Type of supplied Kappa
        enum KMethodType
        {
            mtFluidThermo,
            mtSolidThermo,
            mtDirectionalSolidThermo,
            mtLookup,
            mtFunction
        };


protected:

    // Protected Data

        static const Enum<KMethodType> KMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;

        //- How to get K
        const KMethodType method_;

        //- Name of thermal conductivity field (if looked up from database)
        const word kappaName_;

        //- Name of thermal diffusivity
        const word alphaName_;

        //- Name of the non-isotropic alpha (for directional solidThermo)
        const word alphaAniName_;

        //- Function1 for kappa
        autoPtr<PatchFunction1<scalar>> kappaFunction1_;

        //- Function1 for alpha
        autoPtr<PatchFunction1<scalar>> alphaFunction1_;


public:

    // Constructors

        //- Default construct from patch, using fluidThermo (default)
        //- or specified method
        explicit temperatureCoupledBase
        (
            const fvPatch& patch,
            const KMethodType method = KMethodType::mtFluidThermo
        );

        //- Construct from patch, method type and field names
        temperatureCoupledBase
        (
            const fvPatch& patch,
            const KMethodType method,
            const word& kappaName,
            const word& alphaName,
            const word& alphaAniName
        );

        //- Construct from patch and dictionary
        temperatureCoupledBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );

        //- Construct from patch and temperatureCoupledBase
        temperatureCoupledBase
        (
            const fvPatch& patch,
            const temperatureCoupledBase& base
        );

        //- Copy construct
        temperatureCoupledBase(const temperatureCoupledBase& base);


    //- Destructor
    virtual ~temperatureCoupledBase() = default;


    // Member Functions

        //- Method to obtain K
        word KMethod() const
        {
            return KMethodTypeNames_[method_];
        }

        //- Name of thermal conductivity field
        const word& kappaName() const noexcept
        {
            return kappaName_;
        }

        //- Name of thermal diffusivity field
        const word& alphaName() const noexcept
        {
            return alphaName_;
        }

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        ) = 0;

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<scalar>&,
            const labelList&
        ) = 0;

        //- Given patch temperature calculate corresponding K field
        virtual tmp<scalarField> kappa(const scalarField& Tp) const;

        //- Given patch temperature calculate corresponding alphaEff field
        virtual tmp<scalarField> alpha(const scalarField& Tp) const;

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
