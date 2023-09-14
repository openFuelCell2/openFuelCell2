/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationMode
>
Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationModeNames
({
    { operationMode::fixedPower, "power" },
    { operationMode::fixedHeatFlux, "flux" },
    { operationMode::fixedHeatTransferCoeff, "coefficient" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),  // default method (fluidThermo)
    mode_(fixedHeatFlux),
    Q_(nullptr),
    q_(nullptr),
    h_(nullptr),
    Ta_(nullptr),
    relaxation_(1),
    emissivity_(0),
    qrRelaxation_(1),
    qrName_("undefined-qr"),
    thicknessLayers_(),
    kappaLayers_()
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.get("mode", dict)),
    Q_(nullptr),
    q_(nullptr),
    h_(nullptr),
    Ta_(nullptr),
    relaxation_(dict.getOrDefault<scalar>("relaxation", 1)),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 0)),
    qrRelaxation_(dict.getOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    thicknessLayers_(),
    kappaLayers_()
{
    switch (mode_)
    {
        case fixedPower:
        {
            Q_ = Function1<scalar>::New("Q", dict, &db());
            break;
        }
        case fixedHeatFlux:
        {
            q_ = PatchFunction1<scalar>::New(patch().patch(), "q", dict);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_ = PatchFunction1<scalar>::New(patch().patch(), "h", dict);
            Ta_ = Function1<scalar>::New("Ta", dict, &db());

            if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
            {
                dict.readEntry("kappaLayers", kappaLayers_);

                if (thicknessLayers_.size() != kappaLayers_.size())
                {
                    FatalIOErrorInFunction(dict)
                        << "\n number of layers for thicknessLayers and "
                        << "kappaLayers must be the same"
                        << "\n for patch " << p.name()
                        << " of field " << internalField().name()
                        << " in file " << internalField().objectPath()
                        << exit(FatalIOError);
                }
            }

            break;
        }
    }

    if (qrName_ != "none")
    {
        if (dict.found("qrPrevious"))
        {
            qrPrevious_ = scalarField("qrPrevious", dict, p.size());
        }
        else
        {
            qrPrevious_.resize(p.size(), Zero);
        }
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(rhs, p, iF, mapper),
    temperatureCoupledBase(patch(), rhs),
    mode_(rhs.mode_),
    Q_(rhs.Q_.clone()),
    q_(rhs.q_.clone(patch().patch())),
    h_(rhs.h_.clone(patch().patch())),
    Ta_(rhs.Ta_.clone()),
    relaxation_(rhs.relaxation_),
    emissivity_(rhs.emissivity_),
    qrPrevious_(),
    qrRelaxation_(rhs.qrRelaxation_),
    qrName_(rhs.qrName_),
    thicknessLayers_(rhs.thicknessLayers_),
    kappaLayers_(rhs.kappaLayers_)
{
    if (qrName_ != "none")
    {
        qrPrevious_.resize(mapper.size());
        qrPrevious_.map(rhs.qrPrevious_, mapper);
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& rhs
)
:
    mixedFvPatchScalarField(rhs),
    temperatureCoupledBase(rhs),
    mode_(rhs.mode_),
    Q_(rhs.Q_.clone()),
    q_(rhs.q_.clone(patch().patch())),
    h_(rhs.h_.clone(patch().patch())),
    Ta_(rhs.Ta_.clone()),
    relaxation_(rhs.relaxation_),
    emissivity_(rhs.emissivity_),
    qrPrevious_(rhs.qrPrevious_),
    qrRelaxation_(rhs.qrRelaxation_),
    qrName_(rhs.qrName_),
    thicknessLayers_(rhs.thicknessLayers_),
    kappaLayers_(rhs.kappaLayers_)
{}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(rhs, iF),
    temperatureCoupledBase(patch(), rhs),
    mode_(rhs.mode_),
    Q_(rhs.Q_.clone()),
    q_(rhs.q_.clone(patch().patch())),
    h_(rhs.h_.clone(patch().patch())),
    Ta_(rhs.Ta_.clone()),
    relaxation_(rhs.relaxation_),
    emissivity_(rhs.emissivity_),
    qrPrevious_(rhs.qrPrevious_),
    qrRelaxation_(rhs.qrRelaxation_),
    qrName_(rhs.qrName_),
    thicknessLayers_(rhs.thicknessLayers_),
    kappaLayers_(rhs.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);

    if (q_)
    {
        q_->autoMap(mapper);
    }
    if (h_)
    {
        h_->autoMap(mapper);
    }

    if (qrName_ != "none")
    {
        qrPrevious_.autoMap(mapper);
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const auto& rhs =
        refCast<const externalWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    temperatureCoupledBase::rmap(rhs, addr);


    if (q_)
    {
        q_->rmap(rhs.q_(), addr);
    }
    if (h_)
    {
        h_->rmap(rhs.h_(), addr);
    }

    if (qrName_ != "none")
    {
        qrPrevious_.rmap(rhs.qrPrevious_, addr);
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp(*this);

    const scalarField valueFraction0(valueFraction());
    const scalarField refValue0(refValue());

    scalarField qr(Tp.size(), Zero);
    if (qrName_ != "none")
    {
        qr = lerp
        (
            qrPrevious_,
            patch().lookupPatchField<volScalarField>(qrName_),
            qrRelaxation_
        );
        qrPrevious_ = qr;
    }

    switch (mode_)
    {
        case fixedPower:
        {
            const scalar heatPower =
                Q_->value(this->db().time().timeOutputValue());

            refGrad() = (heatPower/gSum(patch().magSf()) + qr)/kappa(Tp);
            refValue() = 0;
            valueFraction() = 0;

            break;
        }
        case fixedHeatFlux:
        {
            tmp<scalarField> heatFlux =
                q_->value(this->db().time().timeOutputValue());

            refGrad() = (heatFlux + qr)/kappa(Tp);
            refValue() = 0;
            valueFraction() = 0;

            break;
        }
        case fixedHeatTransferCoeff:
        {
            tmp<scalarField> thtcCoeff =
            (
                h_->value(this->db().time().timeOutputValue()) + VSMALL
            );
            const auto& htcCoeff = thtcCoeff();

            scalar totalSolidRes = 0;
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    const scalar l = thicknessLayers_[iLayer];
                    if (kappaLayers_[iLayer] > 0)
                    {
                        totalSolidRes += l/kappaLayers_[iLayer];
                    }
                }
            }

            const scalar Ta =
                Ta_->value(this->db().time().timeOutputValue());

            scalarField hrad(Tp.size(), Zero);

            if (emissivity_ > 0)
            {
                const scalar eSig(emissivity_*sigma.value());
                // Evaluate the radiative flux to the environment
                // from the surface temperature ...
                if (totalSolidRes > 0)
                {
                    // ... including the effect of the solid wall thermal
                    // resistance
                    scalarField TpLambda(htcCoeff/(htcCoeff + 1/totalSolidRes));
                    scalarField Ts(TpLambda*Ta + (1 - TpLambda)*Tp);
                    hrad = eSig*((pow3(Ta) + pow3(Ts)) + Ta*Ts*(Ta + Ts));

                    forAll(hrad, i)
                    {
                        scalar hradTmp0 = hrad[i];
                        scalar TaLambda =
                            (htcCoeff[i] + hradTmp0)
                           /(htcCoeff[i] + hradTmp0 + 1/totalSolidRes);

                        scalar TsiNew = TaLambda*Ta + (1 - TaLambda)*Tp[i];
                        scalar Tsi = Ts[i];

                        while (mag(Tsi - TsiNew)/Tsi > 0.01)
                        {
                            Tsi = TsiNew;
                            scalar hradNew
                            (
                                eSig*((pow3(Ta) + pow3(Tsi)) + Ta*Tsi*(Ta + Tsi))
                            );

                            TaLambda =
                                (htcCoeff[i] + hradNew)
                               /(htcCoeff[i] + hradNew + 1/totalSolidRes);

                            TsiNew = TaLambda*Ta + (1 - TaLambda)*Tp[i];
                        };

                        hrad[i] =
                            eSig*((pow3(Ta) + pow3(Tsi)) + Ta*Tsi*(Ta + Tsi));
                    }
                }
                else
                {
                    hrad = eSig*((pow3(Ta) + pow3(Tp)) + Ta*Tp*(Ta + Tp));
                }
            }

            const scalarField hp(1/(1/(htcCoeff + hrad) + totalSolidRes));

            const scalarField hpTa(hp*Ta);

            const scalarField kappaDeltaCoeffs
            (
                this->kappa(Tp)*patch().deltaCoeffs()
            );

            refGrad() = 0;

            forAll(Tp, i)
            {
                if (qr[i] < 0)
                {
                    const scalar hpmqr = hp[i] - qr[i]/Tp[i];

                    refValue()[i] = hpTa[i]/hpmqr;
                    valueFraction()[i] = hpmqr/(hpmqr + kappaDeltaCoeffs[i]);
                }
                else
                {
                    refValue()[i] = (hpTa[i] + qr[i])/hp[i];
                    valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
                }
            }

            break;
        }
    }

    valueFraction() = lerp(valueFraction0, valueFraction(), relaxation_);
    refValue() = lerp(refValue0, refValue(), relaxation_);

    mixedFvPatchScalarField::updateCoeffs();

    DebugInfo
        << patch().boundaryMesh().mesh().name() << ':' << patch().name() << ':'
        << internalField().name() << " :"
        << " heat transfer rate:" << gSum(kappa(Tp)*patch().magSf()*snGrad())
        << " wall temperature "
        << " min:" << gMin(*this)
        << " max:" << gMax(*this)
        << " avg:" << gAverage(*this) << nl;
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    if (Q_)
    {
        Q_->writeData(os);
    }
    if (q_)
    {
        q_->writeData(os);
    }
    if (h_)
    {
        h_->writeData(os);
    }
    if (Ta_)
    {
        Ta_->writeData(os);
    }

    switch (mode_)
    {
        case fixedHeatTransferCoeff:
        {
            if (relaxation_ < 1)
            {
                os.writeEntry("relaxation", relaxation_);
            }

            if (emissivity_ > 0)
            {
                os.writeEntry("emissivity", emissivity_);
            }

            if (thicknessLayers_.size())
            {
                thicknessLayers_.writeEntry("thicknessLayers", os);
                kappaLayers_.writeEntry("kappaLayers", os);
            }

            break;
        }

        default:
            break;
    }

    os.writeEntry("qr", qrName_);

    if (qrName_ != "none")
    {
        os.writeEntry("qrRelaxation", qrRelaxation_);

        qrPrevious_.writeEntry("qrPrevious", os);
    }

    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
