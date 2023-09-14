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

\*---------------------------------------------------------------------------*/

#include "sorptionWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//- Estimate the y* at the intersection of the two sublayers
static scalar calcYStarLam
(
    const scalar kappa,
    const scalar E,
    const scalar Sc,
    const scalar Sct,
    const scalar Pc
)
{
    scalar ypl = 11;

    for (int iter = 0; iter < 10; ++iter)
    {
        // (F:Eq. 5.5)
        ypl = (log(max(E*ypl, scalar(1)))/kappa + Pc)*Sct/Sc;
    }

    return ypl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<scalarField> sorptionWallFunctionFvPatchScalarField::yPlus() const
{
    // Calculate fields of interest
    const label patchi = patch().index();

    const auto& k = db().lookupObject<volScalarField>(kName_);
    tmp<scalarField> tkwc = k.boundaryField()[patchi].patchInternalField();
    const scalarField& kwc = tkwc.cref();

    const auto& nu = db().lookupObject<volScalarField>(nuName_);
    tmp<scalarField> tnuwc = nu.boundaryField()[patchi].patchInternalField();
    const scalarField& nuwc = tnuwc.cref();

    const volScalarField& y = wallDist::New(internalField().mesh()).y();
    tmp<scalarField> tywc = y.boundaryField()[patchi].patchInternalField();
    const scalarField& ywc = tywc.cref();


    // Calculate the empirical constant given by (Jayatilleke, 1966) (FDC:Eq. 6)
    const scalar Pc =
        9.24*(pow(Sc_/Sct_, 0.75) - 1)*(1 + 0.28*exp(-0.007*Sc_/Sct_));
    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();
    const scalar E = wallCoeffs_.E();

    auto tyPlus = tmp<scalarField>::New(patch().size(), Zero);
    auto& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        // (FDC:Eq. 3)
        const scalar yStar = Cmu25*sqrt(kwc[facei])*ywc[facei]/nuwc[facei];

        // (FDC:Eq. 4)
        const scalar yPlusVis = Sc_*yStar;

        // (FDC:Eq. 5)
        const scalar yPlusLog = Sct_*(log(max(E*yStar, 1 + 1e-4))/kappa + Pc);

        switch (blender_)
        {
            case blenderType::EXPONENTIAL:
            {
                // (FDC:Eq. 2)
                const scalar Gamma =
                    0.01*pow4(Sc_*yStar)/(1 + 5*pow3(Sc_)*yStar);
                const scalar invGamma = scalar(1)/max(Gamma, ROOTVSMALL);

                // (FDC:Eq. 1)
                yPlus[facei] = yPlusVis*exp(-Gamma) + yPlusLog*exp(-invGamma);
                break;
            }

            case blenderType::STEPWISE:
            {
                static const scalar yStarLam =
                    calcYStarLam(kappa, E, Sc_, Sct_, Pc);

                // (F:Eq. 5.3)
                if (yStar < yStarLam)
                {
                    yPlus[facei] = yPlusVis;
                }
                else
                {
                    yPlus[facei] = yPlusLog;
                }
                break;
            }

            case blenderType::BINOMIAL:
            {
                yPlus[facei] =
                    pow
                    (
                        pow(yPlusVis, n_) + pow(yPlusLog, n_),
                        scalar(1)/n_
                    );
                break;
            }

            case blenderType::MAX:
            {
                yPlus[facei] = max(yPlusVis, yPlusLog);
                break;
            }

            case blenderType::TANH:
            {
                const scalar phiTanh = tanh(pow4(0.1*yStar));
                const scalar b1 = yPlusVis + yPlusLog;
                const scalar b2 =
                    pow(pow(yPlusVis, 1.2) + pow(yPlusLog, 1.2), 1.0/1.2);

                yPlus[facei] = phiTanh*b1 + (1 - phiTanh)*b2;
                break;
            }
        }
    }

    return tyPlus;
}


tmp<scalarField> sorptionWallFunctionFvPatchScalarField::flux() const
{
    // Calculate fields of interest
    const label patchi = patch().index();

    const auto& k = db().lookupObject<volScalarField>(kName_);
    tmp<scalarField> tkwc = k.boundaryField()[patchi].patchInternalField();

    const volScalarField& y = wallDist::New(internalField().mesh()).y();
    tmp<scalarField> tywc = y.boundaryField()[patchi].patchInternalField();


    // Calculate mass-transfer coefficient field (FDC:Eqs. 8-9)
    tmp<scalarField> ta =
        laminar_
      ? D_/tywc
      : pow025(wallCoeffs_.Cmu())*sqrt(tkwc)/yPlus();


    // Calculate wall-surface concentration
    const auto& Csurf = *this;

    // Calculate wall-adjacent concentration
    const scalar t = db().time().timeOutputValue();
    tmp<scalarField> tkAbs = kAbsPtr_->value(t);
    tmp<scalarField> tCstar = Csurf/tkAbs;

    // Calculate near-wall cell concentration
    const word& fldName = internalField().name();
    const auto& C = db().lookupObject<volScalarField>(fldName);
    tmp<scalarField> tCwc = C.boundaryField()[patchi].patchInternalField();

    // Return adsorption or absorption/permeation flux
    // normalised by the mass-transfer coefficient (FDC:Fig. 1)
    return (tCstar - tCwc)/ta;
}


void sorptionWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
    os.writeEntryIfDifferent<bool>("laminar", false, laminar_);
    os.writeEntry("Sc", Sc_);
    os.writeEntry("Sct", Sct_);
    os.writeEntryIfDifferent<scalar>("D", -1, D_);
    wallCoeffs_.writeEntries(os);
    os.writeEntryIfDifferent<word>("k", "k", kName_);
    os.writeEntryIfDifferent<word>("nu", "nu", nuName_);

    if (kAbsPtr_)
    {
        kAbsPtr_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sorptionWallFunctionFvPatchScalarField::sorptionWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    wallFunctionBlenders(),
    laminar_(false),
    kAbsPtr_(nullptr),
    Sc_(1),
    Sct_(1),
    D_(-1),
    kName_("k"),
    nuName_("nu"),
    wallCoeffs_()
{}


sorptionWallFunctionFvPatchScalarField::sorptionWallFunctionFvPatchScalarField
(
    const sorptionWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    laminar_(ptf.laminar_),
    kAbsPtr_(ptf.kAbsPtr_.clone(patch().patch())),
    Sc_(ptf.Sc_),
    Sct_(ptf.Sct_),
    D_(ptf.D_),
    kName_(ptf.kName_),
    nuName_(ptf.nuName_),
    wallCoeffs_(ptf.wallCoeffs_)
{}


sorptionWallFunctionFvPatchScalarField::sorptionWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),  // Bypass dictionary constructor
    wallFunctionBlenders(dict, blenderType::STEPWISE, scalar(2)),
    laminar_(dict.getOrDefault<bool>("laminar", false)),
    kAbsPtr_(PatchFunction1<scalar>::New(p.patch(), "kAbs", dict)),
    Sc_(dict.getCheck<scalar>("Sc", scalarMinMax::ge(0))),
    Sct_(dict.getCheck<scalar>("Sct", scalarMinMax::ge(0))),
    D_(dict.getOrDefault<scalar>("D", -1)),
    kName_(dict.getOrDefault<word>("k", "k")),
    nuName_(dict.getOrDefault<word>("nu", "nu")),
    wallCoeffs_(dict)
{
    if (laminar_)
    {
        if (D_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Molecular diffusion coefficient cannot be non-positive. "
                << "D = " << D_
                << exit(FatalIOError);
        }
    }

    if (!kAbsPtr_)
    {
        FatalIOErrorInFunction(dict)
            << "Adsorption or absorption coefficient is not set."
            << exit(FatalIOError);
    }

    if (!this->readGradientEntry(dict) || !this->readValueEntry(dict))
    {
        extrapolateInternal();
        gradient() = Zero;
    }
}


sorptionWallFunctionFvPatchScalarField::sorptionWallFunctionFvPatchScalarField
(
    const sorptionWallFunctionFvPatchScalarField& swfpsf
)
:
    fixedGradientFvPatchScalarField(swfpsf),
    wallFunctionBlenders(swfpsf),
    laminar_(swfpsf.laminar_),
    kAbsPtr_(swfpsf.kAbsPtr_.clone(patch().patch())),
    Sc_(swfpsf.Sc_),
    Sct_(swfpsf.Sct_),
    D_(swfpsf.D_),
    kName_(swfpsf.kName_),
    nuName_(swfpsf.nuName_),
    wallCoeffs_(swfpsf.wallCoeffs_)
{}


sorptionWallFunctionFvPatchScalarField::sorptionWallFunctionFvPatchScalarField
(
    const sorptionWallFunctionFvPatchScalarField& swfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(swfpsf, iF),
    wallFunctionBlenders(swfpsf),
    laminar_(swfpsf.laminar_),
    kAbsPtr_(swfpsf.kAbsPtr_.clone(patch().patch())),
    Sc_(swfpsf.Sc_),
    Sct_(swfpsf.Sct_),
    D_(swfpsf.D_),
    kName_(swfpsf.kName_),
    nuName_(swfpsf.nuName_),
    wallCoeffs_(swfpsf.wallCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sorptionWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedGradientFvPatchScalarField::autoMap(mapper);

    if (kAbsPtr_)
    {
        kAbsPtr_->autoMap(mapper);
    }
}


void sorptionWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const auto& swfptf =
        refCast<const sorptionWallFunctionFvPatchScalarField>(ptf);

    if (kAbsPtr_)
    {
        kAbsPtr_->rmap(swfptf.kAbsPtr_(), addr);
    }
}


void sorptionWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    gradient() = flux()/patch().deltaCoeffs();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void sorptionWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchField<scalar>::write(os);

    writeLocalEntries(os);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    sorptionWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
