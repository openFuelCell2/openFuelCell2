/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd
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
    Foam::compressible::alphatJayatillekeWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    This boundary condition provides a thermal wall function for turbulent
    thermal diffusivity (usually\c alphat) based on the Jayatilleke model.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            alphatJayatillekeWallFunction;

        // Optional entries
        Prt             <scalar>;
        kappa           <scalar>;
        E               <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type   | Reqd  | Deflt
      type      | Type name: alphatJayatillekeWallFunction | word | yes  | -
      Prt       | Turbulent Prandtl number       | scalar | no    | 0.85
      kappa     | von Karman constant            | scalar | no    | 0.41
      E         | Wall roughness parameter       | scalar | no    | 9.8
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

SourceFiles
    alphatJayatillekeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressible_alphatJayatillekeWallFunctionFvPatchScalarField_H
#define compressible_alphatJayatillekeWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
      Class alphatJayatillekeWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatJayatillekeWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;


        // Solution parameters

            static scalar maxExp_;
            static scalar tolerance_;
            static label maxIters_;


    // Private Member Functions

        //- Check the type of the patch
        void checkType();

        //- Return the patch y+
        tmp<scalarField> yPlus
        (
            const compressible::turbulenceModel& turbModel
        ) const;

        //- `P' function
        scalar Psmooth(const scalar Prat) const;

        //- Calculate y+ at the edge of the thermal laminar sublayer
        scalar yPlusTherm
        (
            const scalar P,
            const scalar Prat
        ) const;


public:

    //- Runtime type information
    TypeName("compressible::alphatJayatillekeWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given an
        //- alphatJayatillekeWallFunctionFvPatchScalarField
        //- onto a new patch
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatJayatillekeWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphatJayatillekeWallFunctionFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionFvPatchScalarField&,
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
                new alphatJayatillekeWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
