/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::outletMachNumberPressureFvPatchScalarField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition maintains a certain subsonic Mach number at an
    outlet patch by dynamically adjusting the static outlet pressure. It makes
    it possible, for example, to simulate the flow in a preturbine engine
    exhaust manifold, without resolving details of the flow inside the turbine.
    In general, the flow in a choked nozzle can be non-trivial and expensive
    to simulate.

    This formulation is derived from a simple model of the gas flow through
    a nozzle with fixed geometry. The nozzle flow is assumed to be quasi-steady,
    1D, isentropic and compressible.

    This gives the following general relationship between pressure ratio and
    Mach number in any cross section inside the nozzle:

    \f[
        \frac{p_{tot}}{p}=\left[ 1+ \frac{k-1}{2}\;M^2 \right]^{\frac{k}{k-1}}
    \f]

    where the constant ratio of heat capacities is \f$k=c_p/c_v\f$.
    The Mach number in the cross section is \f$M=V/c\f$, where \f$c\f$ is
    the speed of sound and V is the uniform velocity in the streamwise
    direction.

    Overall pressure difference across the nozzle is
    \f[
        r =  pBack/p_{tot}
    \f]

    When \f$k=1.4\f$, the flow in the nozzle throat becomes choked when
    \f$ r<0.5\f$ and non-choked otherwise. This implementation is not applicable
    when \f$ r>=1 \f$ where backflow would occur.

    The nozzle model assumption locks the relationship between nozzle cross
    sectional areas and Mach numbers. For a choked flow it is only the Mach
    number on the outlet patch, \f$M_{outlet}\f$, that needs to be stated in the
    boundary dictionary.

    Care should be taken however to ensure that the entries in the input
    dictionary and the CFD geometry satisfy the following equation
    \f[
        c1\frac{A_{outlet}}{A_1}=\frac{1}{M_{outlet}}\left[\frac{1+\frac{k-1}{2}
        M_{outlet}^2}{1+\frac{k-1}{2}}\right]^{\frac{k+1}{2(k-1)}}
    \f]
    where \f$c1\f$ compensate for non-uniform outlet profiles, \f$A_{outlet}\f$
    is geometrical outlet patch area and \f$A_1\f$ is assumed nozzle throat
    area.

    In the non-choked case the outlet patch Mach number is calculated as
    \f[
        M_{outlet} =
            \frac{A_1}
            {c1\;A_{outlet}}
            \sqrt{\frac{2}{k-1}\left[r^\frac{2}{k}-r^\frac{k+1}{k} \right]}
    \f]

    The accompanying boundary conditions for velocity should be
    pressureInletOutletVelocity.

    Author: Jens Dahl Kunoy

    Reference:
    \verbatim
       Fox, R.W & McDonald, A. T. (1994).
       Introduction to Fluid Mechanics (4ed SI).
       Wiley
    \endverbatim

Usage

    \table
        Property        | Description              | Required | Default value
        choked          | Defines nozzle conditions| Yes          | None
        relax           | underrelaxation of static pressure| Yes | 0
        M               | outlet Mach number        | Yes (choked) | None
        A1              | Nozzle throat area [m2]   | Yes (non-choked) | 0
        pBack           | Pressure downstream of nozzle| Yes (non-choked) | None
        c1              | Correction factor for non-uniform profiles
                        | No (non-choked) | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            outletMachNumberPressure;
        pBack           101325;
        c1              1;
        A1              0.008;
        relax           0.1;
        choked          false;
        value           uniform 200000;
    }
    \endverbatim

SourceFiles
    outletMachNumberPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef outletMachNumberPressureFvPatchScalarField_H
#define outletMachNumberPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class outletMachNumberPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class outletMachNumberPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Mach number
        scalar M_;

        //- Back pressure
        scalar pBack_;

        //- Model constant
        scalar c1_;

        //- Model constant area input
        scalar A1_;

        //- The name of the flux field
        word phiName_;

        //- The name of the rho field
        word rhoName_;

        //- The name of the U field
        word UName_;

        //- Choked
        Switch choked_;

        //- Relaxation factor
        scalar relax_;


public:

    //- Runtime type information
    TypeName("outletMachNumberPressure");


    // Constructors

        //- Construct from patch and internal field
        outletMachNumberPressureFvPatchScalarField
        (
            const fvPatch& p,
            const DimensionedField<scalar, volMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        outletMachNumberPressureFvPatchScalarField
        (
            const fvPatch& p,
            const DimensionedField<scalar, volMesh>& iF,
            const dictionary& dict
        );

        //- Construct by mapping given
        //- outletMachNumberPressureFvPatchScalarField onto a new patch
        outletMachNumberPressureFvPatchScalarField
        (
            const outletMachNumberPressureFvPatchScalarField& ptf,
            const fvPatch& p,
            const DimensionedField<scalar, volMesh>& iF,
            const fvPatchFieldMapper& mapper
        );

        //- Construct as copy
        outletMachNumberPressureFvPatchScalarField
        (
            const outletMachNumberPressureFvPatchScalarField& tppsf
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new outletMachNumberPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletMachNumberPressureFvPatchScalarField
        (
            const outletMachNumberPressureFvPatchScalarField& tppsf,
            const DimensionedField<scalar, volMesh>& iF
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new outletMachNumberPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
