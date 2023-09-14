/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "outletMappedUniformInletHeatAdditionFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    outletPatchName_(),
    phiName_("phi"),
    Q_(0),
    TMin_(0),
    TMax_(5000)
{}


Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_)
{}


Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    outletPatchName_(dict.get<word>("outletPatch")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    Q_(dict.get<scalar>("Q")),
    TMin_(dict.getOrDefault<scalar>("TMin", 0)),
    TMax_(dict.getOrDefault<scalar>("TMax", 5000))
{}



Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_)
{}



Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::outletMappedUniformInletHeatAdditionFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volScalarField& vsf =
    (
        dynamic_cast<const volScalarField&>(this->internalField())
    );

    const fvPatch& fvp = this->patch();

    label outletPatchID =
        fvp.patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find outlet patch " << outletPatchName_
            << abort(FatalError);
    }

    const fvPatch& outletPatch = fvp.boundaryMesh()[outletPatchID];

    const fvPatchField<scalar>& outletPatchField =
        vsf.boundaryField()[outletPatchID];

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField& outletPatchPhi = phi.boundaryField()[outletPatchID];
    scalar sumOutletPatchPhi = gSum(outletPatchPhi);

    if (sumOutletPatchPhi > SMALL)
    {
        const basicThermo& thermo =
            db().lookupObject<basicThermo>(basicThermo::dictName);

        const scalarField& pp = thermo.p().boundaryField()[outletPatchID];
        const scalarField& pT = thermo.T().boundaryField()[outletPatchID];

        scalar averageOutletField =
            gSum(outletPatchPhi*outletPatchField)/sumOutletPatchPhi;

        const scalarField Cpf(thermo.Cp(pp, pT, outletPatchID));

        scalar totalPhiCp = gSum(outletPatchPhi)*gAverage(Cpf);

        operator==(clamp(averageOutletField + Q_/totalPhiCp, TMin_, TMax_));
    }
    else
    {
        scalar averageOutletField =
            gSum(outletPatch.magSf()*outletPatchField)
           /gSum(outletPatch.magSf());

        operator==(averageOutletField);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::outletMappedUniformInletHeatAdditionFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("outletPatch", outletPatchName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("Q", Q_);
    os.writeEntry("TMin", TMin_);
    os.writeEntry("TMax", TMax_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletMappedUniformInletHeatAdditionFvPatchField
    );
}


// ************************************************************************* //
