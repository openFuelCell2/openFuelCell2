/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd
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
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField

Description
    A thermal wall function for simulation of boiling wall.

    This alpha wall function can handle the following regimes:
        single phase
        subcooled nucleate wall boiling
        transitional boiling
        film boiling.

    The wall function uses a partition method to transfer heat either
    to the liquid or vapor phase. At the moment, this function works
    in a wall temperature fixed mode. i.e, there is no consideration
    for the sudden change of heat transfer coefficient (htc) after
    reaching TDBN (deviation from nucleate boiling temperature).

     References:
    \verbatim
        Numerical simulation of immersion quenching process of an engine cylinder head
        Vedanth Srinivasan, Kil-Min Moon, David Greif, De Ming Wang, Myung-hwan Kim
        Applied Mathematical Modelling 34 (2010) 2111-2128
    \endverbatim


    For the single phase non-boiling regime the standard
    JayatillekeWallFunction is used.

    For the sub-cool nucleate boiling regime the following runtime
    selectable submodels are used:
      - nucleation site density
      - bubble departure frequency
      - bubble departure diameter

    Implements a version of the well-known RPI wall boiling model
    (Kurul & Podowski, 1991). The model implementation is similar to the model
    described by Peltola & Pättikangas (2012) but has been extended with the
    wall heat flux partitioning models.

    References:
    \verbatim
        "On the modeling of multidimensional effects in boiling channels"
        Kurul, N., Podowski, M.Z.,
        ANS Proceedings, National Heat Transfer Conference,
        Minneapolis, Minnesota, USA, July 28-31, 1991,
        ISBN: 0-89448-162-1, pp. 30-40
    \endverbatim

    \verbatim
        "Development and validation of a boiling model for OpenFOAM
        multiphase solver"
        Peltola, J., Pättikangas, T.J.H.,
        CFD4NRS-4 Conference Proceedings, paper 59,
        Daejeon, Korea, September 10-12 2012
    \endverbatim


    The transition boiling regime flux (TBF) is modelled following
    a temperature based linear interpolation between the critical heat flux
    (CHF) and the minimum heat flux (MHF) in such a way that when the wall
    temperature is between the range of TDBN and the Leidenfrost temperature
    (TLeiden) a linear interpolation is used between CHF and MHF.

    Thus, the following models are required:
        LeidenfrostModel
        CHFModel
        CHFSubCoolModel
        MHFModel
        TDNBModel
        filmBoilingModel

    The linear interpolation is as follows:

        TBF = CHF*phi + (1 - phi)*MHF

        where phi:

            phi = wp*(Tw - TDNB)/(TLeiden - TDNB),

        where:
            wp model constant
            Tw wall temperature


    The film boiling regime is applied when Tw is larger than TLeiden. In
    this regime the corrlation from the filmBoilingModel is used for
    calculating the cht from the wall.

    The filmBoilingModel is needed in the vapor field in order to calculate
    the heat transfer to the vapor phase in film boiling regime.


Usage
    \table
        Property     | Description             | Required    | Default value
        phaseType    | 'vapor' or 'liquid'     | yes         |
        relax        |wall boiling model relaxation| yes     |
        Prt          | inherited from alphatPhaseChangeJayatillekeWallFunction
        Cmu          | inherited from alphatPhaseChangeJayatillekeWallFunction
        kappa        | inherited from alphatPhaseChangeJayatillekeWallFunction
        E            | inherited from alphatPhaseChangeJayatillekeWallFunction
        dmdt         | phase change mass flux  | no         |
        value        | initial alphat value    | yes         |

        if phaseType 'vapor':

        partitioningModel|                     | yes        |
        filmBoilingModel |                     | yes        |
        LeidenfrostModel |                     | yes        |

        if phaseType 'liquid':

        partitioningModel|                     | yes        |
        nucleationSiteModel|                     | yes      |
        departureDiamModel|                     | yes       |
        departureFreqModel|                     | yes       |
        K                | bubbles area constant| no       | 4

        LeidenfrostModel |                      | no       |
        CHFModel         |                      | no       |
        CHFSubCoolModel  |                      | no       |
        MHFModel         |                      | no       |
        TDNBModel        |                      | no       |
        filmBoilingModel |                      | no       |
        wp               |                      | no       | 1
    \endtable

    NOTE: Runtime selectabale submodels may require model specific entries

    Example usage:
    \verbatim
    hotWall
    {
        type            compressible::alphatWallBoilingWallFunction;
        phaseType       liquid;
        Prt             0.85;
        Cmu             0.09;
        kappa           0.41;
        E               9.8;
        relax           0.1;
        dmdt            uniform 0;
        partitioningModel
        {
            type        Lavieville;
            alphaCrit   0.2;
        }
        nucleationSiteModel
        {
            type        LemmertChawla;
        }
        departureDiamModel
        {
            type        TolubinskiKostanchuk;
        }
        departureFreqModel
        {
            type        Cole;
        }

        LeidenfrostModel
        {
            type        Spiegler;
            Tcrit       647;
        }
        CHFModel
        {
            type        Zuber;
        }
        CHFSubCoolModel
        {
            type        HuaXu;
            Kburn       0.5;
        }
        MHFModel
        {
            type        Jeschar;
            Kmhf        1;
        }
        TDNBModel
        {
            type        Schroeder;
        }
        filmBoilingModel
        {
            type        Bromley;
        }
        value           uniform 0.01;
    \endverbatim

See also
    Foam::alphatPhaseChangeJayatillekeWallFunctionFvPatchField

SourceFiles
    alphatWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressible_alphatWallBoilingWallFunctionFvPatchScalarField_H
#define compressible_alphatWallBoilingWallFunctionFvPatchScalarField_H

#include "Function1.H"

#include "alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H"
#include "partitioningModel.H"
#include "nucleationSiteModel.H"
#include "departureDiameterModel.H"
#include "departureFrequencyModel.H"

#include "LeidenfrostModel.H"
#include "filmBoilingModel.H"
#include "CHFModel.H"
#include "CHFSubCoolModel.H"
#include "MHFModel.H"
#include "TDNBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class alphatWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatWallBoilingWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
{
public:

    // Data types

        //- Enumeration listing the possible operational modes
        enum phaseType
        {
            vaporPhase,
            liquidPhase
        };


private:

    // Private data

        //- Enumeration of regimes per face
        enum regimeType
        {
            subcool,
            transient,
            film,
            nonBoiling
        };

        //- name of the other phase (vapor/liquid phase)
        word otherPhaseName_;

        //- Heat source type names
        static const Enum<phaseType> phaseTypeNames_;

        //- Heat source type
        phaseType phaseType_;

        //- dmdt relaxationFactor
        autoPtr<Function1<scalar>> relax_;

        //- Patch face area by cell volume
        scalarField AbyV_;

        // Sub-cooling nucleating boiling

            //- Convective turbulent thermal diffusivity
            scalarField alphatConv_;

            //- Departure diameter field
            scalarField dDep_;

            //- Quenching surface heat flux
            scalarField qq_;

            //- Model constant for area of bubbles
            scalar K_;

            //- Run-time selected heat flux partitioning model
            autoPtr<wallBoilingModels::partitioningModel>
                partitioningModel_;

            //- Run-time selected nucleation site density model
            autoPtr<wallBoilingModels::nucleationSiteModel>
                nucleationSiteModel_;

            //- Run-time selected bubble departure diameter model
            autoPtr<wallBoilingModels::departureDiameterModel>
                departureDiamModel_;

            //- Run-time selected bubble departure frequency model
            autoPtr<wallBoilingModels::departureFrequencyModel>
                departureFreqModel_;


        // Film boiling model

            //- Run-time selected for filmBoiling model
            autoPtr<wallBoilingModels::filmBoilingModel>
                filmBoilingModel_;

        // Transition boiling model

            //- Run-time selected for Leidenfrost temperature
            autoPtr<wallBoilingModels::LeidenfrostModel>
                LeidenfrostModel_;

            //- Run-time selected for CHF
            autoPtr<wallBoilingModels::CHFModel> CHFModel_;

            //- Run-time selected for CHF sub-cool
            autoPtr<wallBoilingModels::CHFSubCoolModel> CHFSoobModel_;

            //- Run-time selected for MHF
            autoPtr<wallBoilingModels::MHFModel> MHFModel_;

             //- Run-time selected for MHF
            autoPtr<wallBoilingModels::TDNBModel> TDNBModel_;

            //- Wetting parameter for transient boiling
            scalar wp_;


public:

    //- Runtime type information
    TypeName("compressible::alphatWallBoilingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatWallBoilingWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatWallBoilingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
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
                new alphatWallBoilingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        using alphatPhaseChangeWallFunctionFvPatchScalarField::dmdt;

        //- Is there phase change mass transfer for this phasePair
        virtual bool activePhasePair(const phasePairKey&) const;

        //- Return the rate of phase-change for specific phase pair
        virtual const scalarField& dmdt(const phasePairKey&) const;

        //- Return the rate of phase-change for specific phase pair
        virtual const scalarField& mDotL(const phasePairKey&) const;

        //- Return the departure diameter field
        const scalarField& dDeparture() const
        {
            return dDep_;
        }

        //- Return the quenching surface heat flux [W/m2]
        const scalarField& qq() const
        {
            return qq_;
        }

        //- Return the evaporation surface heat flux [W/m2]
        tmp<scalarField> qe() const
        {
            return mDotL_/AbyV_;
        }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
