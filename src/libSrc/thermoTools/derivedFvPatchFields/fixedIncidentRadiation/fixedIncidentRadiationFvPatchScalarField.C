/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "fixedIncidentRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "constants.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),  // default method (fluidThermo)
    qrIncident_(p.size(), Zero)
{}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    qrIncident_(psf.qrIncident_)
{}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),  // Bypass dictionary constructor
    temperatureCoupledBase(patch(), dict),
    qrIncident_("qrIncident", dict, p.size())
{
    if (!this->readGradientEntry(dict) || !this->readValueEntry(dict))
    {
        extrapolateInternal();
        gradient() = Zero;
    }
}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    qrIncident_(psf.qrIncident_)
{}


Foam::radiation::fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    temperatureCoupledBase(patch(), ptf),
    qrIncident_(ptf.qrIncident_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    temperatureCoupledBase::autoMap(m);
    qrIncident_.autoMap(m);
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(psf, addr);

    const fixedIncidentRadiationFvPatchScalarField& thftpsf =
        refCast<const fixedIncidentRadiationFvPatchScalarField>
        (
            psf
        );

    temperatureCoupledBase::rmap(thftpsf, addr);
    qrIncident_.rmap(thftpsf.qrIncident_, addr);
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField intFld(patchInternalField());

    const radiation::radiationModel& radiation =
        db().lookupObject<radiation::radiationModel>("radiationProperties");

    scalarField emissivity
    (
        radiation.absorptionEmission().e()().boundaryField()[patch().index()]
    );

    gradient() =
        emissivity
       *(
            qrIncident_
          - physicoChemical::sigma.value()*pow4(*this)
        )/kappa(*this);

    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar qr = gSum(kappa(*this)*gradient()*patch().magSf());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " -> "
            << " radiativeFlux:" << qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::radiation::fixedIncidentRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    temperatureCoupledBase::write(os);
    qrIncident_.writeEntry("qrIncident", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedIncidentRadiationFvPatchScalarField
    );
}
}

// ************************************************************************* //
