/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::compressible::turbulentTemperatureRadCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer,
    suitable for multiregion cases. The optional specification of
    thin thermal layer resistances can be made using the entries
    \c thicknessLayers and \c kappaLayers.

    The thermal conductivity \c kappa can be obtained from a variety
    of sources, as explained in the \c temperatureCoupledBase class.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type                compressible::turbulentTemperatureRadCoupledMixed;

        // Optional entries
        Tnbr                <word>;
        qrNbr               <word>;
        qr                  <word>;
        logInterval         <scalar>;
        thermalInertia      <bool>;
        verbose             <bool>;
        prefix              <word>;

        // Conditional entries

            // Option-1
            thicknessLayers <scalarList>;
            kappaLayers     <scalarList>;

            // Option-2
            thicknessLayer  <PatchFunction1<scalar>>;
            kappaLayer      <PatchFunction1<scalar>>;

        // Inherited entries
        kappaMethod         <word>;
        kappa               <word>;
        ...
    }
    \endverbatim

    For example:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureRadCoupledMixed;
        Tnbr            T;
        qrNbr           qr; // or none. Name of qr field on neighbour region
        qr              qr; // or none. Name of qr field on this region
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        thermalInertia  false/true;
        kappaMethod     lookup;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type  | Reqd | Deflt
      type   | Type name: compressible::turbulentTemperatureRadCoupledMixed <!--
               -->                                      | word  | yes  | -
      Tnbr     | Name of the nbr temperature field      | word  | no   | T
      qrNbr    | Name of the radiative flux in the nbr region | word | no | none
      qr       | Name of the radiative flux in this region | word | no | none
      logInterval | Log-file output frequency [s]          | scalar | no  | -1
      thermalInertia | Flag to add thermal inertia to wall node | bool | no <!--
               -->                                                     | false
      verbose | Flag to enable verbose output with extra fields | bool | no <!--
              -->                                                      | false
      prefix  | Name of output-field prefix (in verbose mode)         <!--
              -->                                       | word | no | multiWorld
      thicknessLayers | List of thicknesses per layer [m] | scalarList <!--
              -->                                             | choice | -
      kappaLayers  | List of thermal conductivites per layer [W/m/K]   <!--
              -->                                | scalarList | choice | -
      thicknessLayer | Single thickness of layer [m]                   <!--
              -->                   | PatchFunction1\<scalar\>| choice | -
      kappaLayer     | Corresponding thermal conductivity [W/m/K]      <!--
              -->                   | PatchFunction1\<scalar\>| choice | -
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchFields.H \endlink
      - \link temperatureCoupledBase.H \endlink
      - \link mappedPatchFieldBase.H \endlink
      - \link writeFile.H \endlink

Note
  - Needs to be on underlying \c mapped(Wall)FvPatch.

SourceFiles
    turbulentTemperatureRadCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_turbulentTemperatureRadCoupledMixedFvPatchScalarField_H
#define Foam_turbulentTemperatureRadCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "scalarList.H"
#include "mappedPatchFieldBase.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class turbulentTemperatureRadCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureRadCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase,
    public mappedPatchFieldBase<scalar>,
    public functionObjects::writeFile
{
    // Private Data

        //- Name of temperature field on the neighbour region
        const word TnbrName_;

        //- Name of the radiative heat flux in the neighbour region
        const word qrNbrName_;

        //- Name of the radiative heat flux in this region
        const word qrName_;

        //- Thickness of layers (either scalarList or a single PatchFunction1)
        scalarList thicknessLayers_;
        autoPtr<PatchFunction1<scalar>> thicknessLayer_;

        //- Conductivity of layers
        scalarList kappaLayers_;
        autoPtr<PatchFunction1<scalar>> kappaLayer_;

        //- Log-file output frequency
        scalar logInterval_;

        //- Index of previous execution
        label executionIndex_;

        //- Flag to enable thermal inertia term
        const Switch thermalInertia_;

        //- Flag to enable verbose output with extra fields
        const bool verbose_;

        //- Name of output-field prefix (in verbose mode)
        const word prefix_;


     // Private Functions

        //- Return local alphaSfDelta
        tmp<scalarField> alphaSfDelta() const;

        //- Return delta enthalpy between regions
        tmp<scalarField> deltaH() const;

        //- Return the sum of deltaCoeff*alpha from nbr and local
        tmp<scalarField> beta() const;

        //- Calculate coefficients for assembly matrix
        tmp<Field<scalar>> coeffs
        (
            fvMatrix<scalar>& matrix,
            const Field<scalar>&,
            const label
        ) const;

        //- Return requested field from the object registry
        //- or create+register the field to the object registry
        volScalarField& getOrCreateField(const word& fieldName) const;

        //- Store heat-transfer coeff fields to facilitate postprocessing
        void storeHTCFields
        (
            const word& prefix,
            const scalarField& htc,
            const scalarField& htcPatch
        ) const;

        //- Flag to allow writing to file
        bool writeFile();


        // functionObjects::writeFile

            //- Write file header information
            void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureRadCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //- new patch
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const
            turbulentTemperatureRadCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const turbulentTemperatureRadCoupledMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureRadCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const turbulentTemperatureRadCoupledMixedFvPatchScalarField&,
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
                new turbulentTemperatureRadCoupledMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

    // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<scalar>&,
            const labelList&
        );


    // Evaluation

        //- Given patch temperature calculate corresponding K field. Override
        //- temperatureCoupledBase::kappa to includes effect of any
        //- explicit kappaThickness
        virtual tmp<scalarField> kappa(const scalarField& Tp) const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Manipulate matrix
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>& m,
            const label iMatrix,
            const direction cmpt
        );


    // I-O

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
