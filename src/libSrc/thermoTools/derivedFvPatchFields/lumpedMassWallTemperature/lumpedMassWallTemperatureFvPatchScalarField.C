/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "lumpedMassWallTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),  // default method (fluidThermo)
    Cp_(0.0),
    mass_(0.0),
    curTimeIndex_(-1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    Cp_(ptf.Cp_),
    mass_(ptf.mass_),
    curTimeIndex_(-1)
{}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    Cp_(dict.get<scalar>("Cp")),
    mass_(dict.get<scalar>("mass")),
    curTimeIndex_(-1)
{
    fvPatchFieldBase::readDict(dict);
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedMassWallTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
}


void Foam::lumpedMassWallTemperatureFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const lumpedMassWallTemperatureFvPatchScalarField& tiptf =
        refCast
        <
            const lumpedMassWallTemperatureFvPatchScalarField
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
}


void Foam::lumpedMassWallTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }

    // Calculate heat flux in or out the wall
    scalarField& Tp(*this);
    const scalarField& magSf = patch().magSf();

    const scalar deltaT(db().time().deltaTValue());

    tmp<scalarField> tkappa(kappa(Tp));

    const scalarField q(tkappa.ref()*snGrad());

    // Total heat in or out of the wall
    const scalar Q = gSum(q*magSf);

    Tp += -(Q/mass_/Cp_)*deltaT;

    refGrad() = 0.0;
    refValue() = Tp;
    valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qin(0);
        scalar Qout(0);

        forAll(q, facei)
        {
            if (q[facei] > 0.0) // out the wall
            {
                Qout += q[facei]*magSf[facei];
            }
            else if (q[facei] < 0.0) // into the wall
            {
                Qin += q[facei]*magSf[facei];
            }
        }

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << " Qin [W]:" << Qin
            << " Qout [W]:" << Qout
            << endl;
    }

    curTimeIndex_ = this->db().time().timeIndex();
}


void Foam::lumpedMassWallTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    temperatureCoupledBase::write(os);

    os.writeEntry("Cp", Cp_);
    os.writeEntry("mass", mass_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        lumpedMassWallTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
