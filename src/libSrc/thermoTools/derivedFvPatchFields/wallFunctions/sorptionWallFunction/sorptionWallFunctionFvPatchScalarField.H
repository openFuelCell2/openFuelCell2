/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::sorptionWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    The \c sorptionWallFunction is a wall boundary condition to
    specify scalar/concentration gradient for turbulent and laminar flows.

    The governing equation of the boundary condition is:

    \f[
        \nabla C = \frac{C^* - C_p}{\Delta_y} = \frac{F}{a \Delta_y}
    \f]

    with

    \f[
        C^* = \frac{C_{surf}}{K}
    \f]

    and with the mass-transfer coefficient is calculated for turbulent flows

    \f[
        a = \frac{C_\mu^{0.25} k_p^{0.5}}{y^+_{blended}}
    \f]

    or for laminar-flow and molecular-diffusion-only states

    \f[
        a = \frac{D_m}{y_1}
    \f]

    where
    \vartable
      \nabla C      | Gradient of concentration
      C^*           | Wall-adjacent concentration
      C_p           | Near-wall cell concentration
      \Delta_y      | First-cell centre wall distance
      F             | Flux of concentration
      a             | Mass-transfer coefficient
      C_{surf}      | Wall-surface concentration
      K             | Adsorption or absorption/permeation coefficient
      C_\mu         | Empirical model coefficient
      k_p           | Turbulent kinetic energy in near-wall cell
      y^+_{blended} | Non-dimensional blended near-wall cell height
      D_m           | Molecular-diffusion coefficient
      y_1           | First-cell centre wall distance
    \endvartable

    Required fields:
    \verbatim
      x    | Arbitrary scalar field, e.g. species, passive scalars etc.
    \endverbatim

    Reference:
    \verbatim
        Standard model for exponential blending (tag:FDC):
            Foat, T., Drodge, J., Charleson, A., Whatmore, B.,
            Pownall, S., Glover, P., ... & Marr, A. (2022).
            Predicting vapour transport from semi-volatile organic
            compounds concealed within permeable packaging.
            International Journal of Heat and Mass Transfer, 183, 122012.
            DOI:10.1016/j.ijheatmasstransfer.2021.122012

        Standard model for stepwise blending (tag:F):
            Foat, T. (2021).
            Modelling vapour transport in indoor environments for
            improved detection of explosives using dogs.
            Doctoral dissertation. University of Southampton.
            URI:http://eprints.soton.ac.uk/id/eprint/456709
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            sorptionWallFunction;
        Sc              <scalar>;
        Sct             <scalar>;
        kAbs            <PatchFunction1<scalar>>;

        // Optional entries
        laminar         <bool>;
        D               <scalar>;
        kName           <word>;
        nuName          <word>;

        // Inherited entries
        Cmu             <scalar>;
        kappa           <scalar>;
        E               <scalar>;
        blending        <word>;
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type | Reqd | Deflt
      type     | Type name: sorptionWallFunction  | word | yes  | -
      Sc       | Schmidt number                 | scalar | yes  | -
      Sct      | Turbulent Schmidt number       | scalar | yes  | -
      kAbs | Adsorption or absorption/permeation coefficient <!--
               -->            | PatchFunction1\<scalar\> | yes  | -
      laminar  | Flag to calculate mass-transfer coefficient under the <!--
               --> laminar-flow or molecular-diffusion-only states <!--
               -->                                | bool | no   | false
      kName    | Name of operand turbulent kinetic energy field | word | no | k
      nuName   | Name of operand kinematic viscosity field      | word | no | nu
    \endtable

    The inherited entries are elaborated in:
      - \link fixedGradientFvPatchField.H \endlink
      - \link wallFunctionCoefficients.H \endlink
      - \link wallFunctionBlenders.H \endlink
      - \link PatchFunction1.H \endlink

SourceFiles
    sorptionWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sorptionWallFunctionFvPatchScalarFields_H
#define Foam_sorptionWallFunctionFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "wallFunctionCoefficients.H"
#include "wallFunctionBlenders.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class sorptionWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sorptionWallFunctionFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    private wallFunctionBlenders
{
        //- Flag to calculate mass-transfer coefficient
        //- under the laminar-flow or molecular-diffusion-only states
        bool laminar_;

        //- Adsorption or absorption/permeation coefficient
        autoPtr<PatchFunction1<scalar>> kAbsPtr_;

        //- Schmidt number
        scalar Sc_;

        //- Turbulent Schmidt number
        scalar Sct_;

        //- Molecular diffusion coefficient
        scalar D_;

        //- Name of operand turbulent kinetic energy field
        word kName_;

        //- Name of operand kinematic viscosity field
        word nuName_;

        //- Standard wall-function coefficients
        wallFunctionCoefficients wallCoeffs_;


    // Private Member Functions

        //- Return the non-dimensional near-wall cell height field
        //- blended between the viscous and inertial sublayers
        tmp<scalarField> yPlus() const;

        //- Return the flux normalised by the mass-transfer coefficient
        tmp<scalarField> flux() const;

        //- Write local wall-function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("sorptionWallFunction");


    // Constructors

        //- Construct from patch and internal field
        sorptionWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sorptionWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- sorptionWallFunctionFvPatchScalarField onto
        //- a new patch
        sorptionWallFunctionFvPatchScalarField
        (
            const sorptionWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sorptionWallFunctionFvPatchScalarField
        (
            const sorptionWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sorptionWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sorptionWallFunctionFvPatchScalarField
        (
            const sorptionWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new sorptionWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

    // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        );


    // Evaluation

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


    // I-O

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
