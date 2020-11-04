/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "alphatWallBoilingWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "saturationModel.H"
#include "wallFvPatch.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::compressible::
    alphatWallBoilingWallFunctionFvPatchScalarField::phaseType
>
Foam::compressible::
alphatWallBoilingWallFunctionFvPatchScalarField::phaseTypeNames_
{
    { phaseType::vaporPhase, "vapor" },
    { phaseType::liquidPhase, "liquid" },
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF),
    otherPhaseName_("vapor"),
    phaseType_(liquidPhase),
    relax_(),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    qq_(p.size(), 0),
    K_(4),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr),
    filmBoilingModel_(nullptr),
    LeidenfrostModel_(nullptr),
    CHFModel_(nullptr),
    CHFSoobModel_(nullptr),
    MHFModel_(nullptr),
    TDNBModel_(nullptr),
    wp_(1)
{
    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    otherPhaseName_(dict.get<word>("otherPhase")),
    phaseType_(phaseTypeNames_.get("phaseType", dict)),
    relax_(Function1<scalar>::New("relax", dict)),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    qq_(p.size(), 0),
    K_(4),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr),
    filmBoilingModel_(nullptr),
    LeidenfrostModel_(nullptr),
    CHFModel_(nullptr),
    CHFSoobModel_(nullptr),
    MHFModel_(nullptr),
    TDNBModel_(nullptr),
    wp_(1)
{

    // Check that otherPhaseName != this phase
    if (internalField().group() == otherPhaseName_)
    {
        FatalErrorInFunction
            << "otherPhase should be the name of the vapor phase that "
            << "corresponds to the liquid base of vice versa" << nl
            << "This phase: " << internalField().group() << nl
            << "otherPhase: " << otherPhaseName_
            << abort(FatalError);
    }

    switch (phaseType_)
    {
        case vaporPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            dmdt_ = 0;

            break;
        }
        case liquidPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            nucleationSiteModel_ =
                wallBoilingModels::nucleationSiteModel::New
                (
                    dict.subDict("nucleationSiteModel")
                );

            departureDiamModel_ =
                wallBoilingModels::departureDiameterModel::New
                (
                    dict.subDict("departureDiamModel")
                );

            departureFreqModel_ =
                wallBoilingModels::departureFrequencyModel::New
                (
                    dict.subDict("departureFreqModel")
                );

            const dictionary* LeidenfrostDict =
                dict.findDict("LeidenfrostModel");

            if (LeidenfrostDict)
            {
                LeidenfrostModel_ =
                    wallBoilingModels::LeidenfrostModel::New(*LeidenfrostDict);
            }

            const dictionary* CHFDict = dict.findDict("CHFModel");

            if (CHFDict)
            {
                CHFModel_ =
                    wallBoilingModels::CHFModel::New(*CHFDict);
            }

            const dictionary* HFSubCoolDict = dict.findDict("CHFSubCoolModel");

            if (HFSubCoolDict)
            {
                CHFSoobModel_ =
                    wallBoilingModels::CHFSubCoolModel::New(*HFSubCoolDict);
            }

            const dictionary* MHFDict = dict.findDict("MHFModel");

            if (MHFDict)
            {
                MHFModel_ =
                    wallBoilingModels::MHFModel::New(*MHFDict);
            }

            const dictionary* TDNBDict = dict.findDict("TDNBModel");

            if (TDNBDict)
            {
                TDNBModel_ =
                    wallBoilingModels::TDNBModel::New(*TDNBDict);
            }

            const dictionary* filmDict = dict.findDict("filmBoilingModel");

            if (filmDict)
            {
                filmBoilingModel_ =
                    wallBoilingModels::filmBoilingModel::New(*filmDict);
            }

            if (dict.found("dDep"))
            {
                dDep_ = scalarField("dDep", dict, p.size());
            }

            dict.readIfPresent("K", K_);

            dict.readIfPresent("wp", wp_);

            if (dict.found("qQuenching"))
            {
                qq_ = scalarField("qQuenching", dict, p.size());
            }

            break;
        }
    }

    if (dict.found("alphatConv"))
    {
        alphatConv_ = scalarField("alphatConv", dict, p.size());
    }

    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_.clone()),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_, mapper),
    dDep_(psf.dDep_, mapper),
    qq_(psf.qq_, mapper),
    K_(psf.K_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    filmBoilingModel_(psf.filmBoilingModel_),
    LeidenfrostModel_(psf.LeidenfrostModel_),
    CHFModel_(psf.CHFModel_),
    CHFSoobModel_(psf.CHFSoobModel_),
    MHFModel_(psf.MHFModel_),
    TDNBModel_(psf.TDNBModel_),
    wp_(psf.wp_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_.clone()),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    qq_(psf.qq_),
    K_(psf.K_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    filmBoilingModel_(psf.filmBoilingModel_),
    LeidenfrostModel_(psf.LeidenfrostModel_),
    CHFModel_(psf.CHFModel_),
    CHFSoobModel_(psf.CHFSoobModel_),
    MHFModel_(psf.MHFModel_),
    TDNBModel_(psf.TDNBModel_),
    wp_(psf.wp_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_.clone()),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    qq_(psf.qq_),
    K_(psf.K_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    filmBoilingModel_(psf.filmBoilingModel_),
    LeidenfrostModel_(psf.LeidenfrostModel_),
    CHFModel_(psf.CHFModel_),
    CHFSoobModel_(psf.CHFSoobModel_),
    MHFModel_(psf.MHFModel_),
    TDNBModel_(psf.TDNBModel_),
    wp_(psf.wp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool alphatWallBoilingWallFunctionFvPatchScalarField::
activePhasePair(const phasePairKey& phasePair) const
{
    if (phasePair == phasePairKey(otherPhaseName_, internalField().group()))
    {
        return true;
    }
    else
    {
        return false;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarField::
dmdt(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return dmdt_;
    }
    else
    {
        FatalErrorInFunction
            << " dmdt requested for invalid phasePair!"
            << abort(FatalError);

        return dmdt_;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarField::
mDotL(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return mDotL_;
    }
    else
    {
        FatalErrorInFunction
            << " mDotL requested for invalid phasePair!"
            << abort(FatalError);

        return mDotL_;
    }
}

void alphatWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    // Check that partitioningModel has been constructed
    if (!partitioningModel_)
    {
        FatalErrorInFunction
            << "partitioningModel has not been constructed!"
            << abort(FatalError);
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const saturationModel& satModel =
        db().lookupObject<saturationModel>("saturationModel");

    const label patchi = patch().index();

    const scalar t = this->db().time().timeOutputValue();
    scalar relax = relax_->value(t);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            const phaseModel& vapor
            (
                fluid.phases()[internalField().group()]
            );

            const fvPatchScalarField& hewv =
                vapor.thermo().he().boundaryField()[patchi];

            // Vapor Liquid phase fraction at the wall
            const scalarField vaporw(vapor.boundaryField()[patchi]);

            // NOTE! Assumes 1-thisPhase for liquid fraction in
            // multiphase simulations
            const scalarField fLiquid
            (
                partitioningModel_->fLiquid(1-vaporw)
            );

            const tmp<scalarField> talphaw = vapor.thermo().alpha(patchi);
            const scalarField& alphaw = talphaw();

            const scalarField heSnGrad(max(hewv.snGrad(), scalar(1e-16)));

            // Convective thermal diffusivity for single phase
            const scalarField alphatv(calcAlphat(*this));

            forAll(*this, i)
            {
                this->operator[](i) =
                (
                    (1 - fLiquid[i])*(alphatv[i])
                   /max(vaporw[i], scalar(1e-8))
                );
            }

            if (debug)
            {
                Info<< "alphat for vapour : " << nl << endl;

                Info<< "  alphatEffv: " << gMin(vaporw*(*this + alphaw))
                    << " - " << gMax(vaporw*(*this + alphaw)) << endl;

                const scalarField qEff(vaporw*(*this + alphaw)*hewv.snGrad());

                scalar Qeff = gSum(qEff*patch().magSf());
                Info<< " Effective heat transfer rate to vapor:" << Qeff
                    << nl << endl;
            }
            break;
        }
        case liquidPhase:
        {
            // Check that nucleationSiteModel has been constructed
            if (!nucleationSiteModel_)
            {
                FatalErrorInFunction
                    << "nucleationSiteModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that departureDiameterModel has been constructed
            if (!departureDiamModel_)
            {
                FatalErrorInFunction
                    << "departureDiameterModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that nucleationSiteModel has been constructed
            if (!departureFreqModel_)
            {
                FatalErrorInFunction
                    << "departureFrequencyModel has not been constructed!"
                    << abort(FatalError);
            }

            const phaseModel& liquid
            (
                fluid.phases()[internalField().group()]
            );

            const phaseModel& vapor(fluid.phases()[otherPhaseName_]);

            // Retrieve turbulence properties from models
            const phaseCompressibleTurbulenceModel& turbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                );
            const phaseCompressibleTurbulenceModel& vaporTurbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        vapor.name()
                    )
                );

            const tmp<scalarField> tnutw = turbModel.nut(patchi);

            const scalar Cmu25(pow025(Cmu_));

            const scalarField& y = turbModel.y()[patchi];

            const tmp<scalarField> tmuw = turbModel.mu(patchi);
            const scalarField& muw = tmuw();

            const tmp<scalarField> talphaw = liquid.thermo().alphahe(patchi);
            const scalarField& alphaw = talphaw();

            const tmp<volScalarField> tk = turbModel.k();
            const volScalarField& k = tk();
            const fvPatchScalarField& kw = k.boundaryField()[patchi];

            const fvPatchVectorField& Uw =
                turbModel.U().boundaryField()[patchi];
            const scalarField magUp(mag(Uw.patchInternalField() - Uw));
            const scalarField magGradUw(mag(Uw.snGrad()));

            const fvPatchScalarField& rhow =
                turbModel.rho().boundaryField()[patchi];


            const fvPatchScalarField& Tw =
                liquid.thermo().T().boundaryField()[patchi];
            const scalarField Tc(Tw.patchInternalField());

            const scalarField uTau(Cmu25*sqrt(kw));

            const scalarField yPlus(uTau*y/(muw/rhow));

            const scalarField Pr(muw/alphaw);

            // Molecular-to-turbulent Prandtl number ratio
            const scalarField Prat(Pr/Prt_);

            // Thermal sublayer thickness
            const scalarField P(this->Psmooth(Prat));

            const scalarField yPlusTherm(this->yPlusTherm(P, Prat));

            const fvPatchScalarField& rhoVaporw =
                vaporTurbModel.rho().boundaryField()[patchi];

            tmp<volScalarField> tCp = liquid.thermo().Cp();
            const volScalarField& Cp = tCp();
            const fvPatchScalarField& Cpw = Cp.boundaryField()[patchi];

            // Saturation temperature
            const tmp<volScalarField> tTsat =
                satModel.Tsat(liquid.thermo().p());

            const volScalarField& Tsat = tTsat();
            const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);
            const scalarField Tsatc(Tsatw.patchInternalField());

            const fvPatchScalarField& pw =
                liquid.thermo().p().boundaryField()[patchi];

            const fvPatchScalarField& hew =
                liquid.thermo().he().boundaryField()[patchi];

            const scalarField hw
            (
                liquid.thermo().he().member() == "e"
              ? hew.patchInternalField() + pw/rhow.patchInternalField()
              : hew.patchInternalField()
            );

            const scalarField L
            (
                vapor.thermo().he().member() == "e"
              ? vapor.thermo().he(pw, Tsatc, patchi) + pw/rhoVaporw - hw
              : vapor.thermo().he(pw, Tsatc, patchi) - hw
            );

            // Liquid phase fraction at the wall
            const scalarField liquidw(liquid.boundaryField()[patchi]);

            const scalarField fLiquid(partitioningModel_->fLiquid(liquidw));

            // Liquid temperature at y+=250 is estimated from logarithmic
            // thermal wall function (Koncar, Krepper & Egorov, 2005)
            const scalarField Tplus_y250
            (
                Prt_*(log(E_*250)/kappa_ + P)
            );
            const scalarField Tplus(Prt_*(log(E_*yPlus)/kappa_ + P));
            scalarField Tl(Tw - (Tplus_y250/Tplus)*(Tw - Tc));
            Tl = max(Tc - 40, Tl);

            // Film, transient boiling regimes
            scalarField Qtb(this->size(), 0);
            scalarField tDNB(this->size(), GREAT);
            scalarField TLeiden(this->size(), GREAT);
            scalarField htcFilmBoiling(this->size(), 0);

            if
            (
                CHFModel_.valid()
                && CHFSoobModel_.valid()
                && TDNBModel_.valid()
                && MHFModel_.valid()
                && LeidenfrostModel_.valid()
                && filmBoilingModel_.valid()
            )
            {

                const scalarField CHF
                (
                    CHFModel_->CHF
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    )
                );

                // Effect of sub-cooling to the CHF in saturated conditions
                const scalarField CHFSubCool
                (
                    CHFSoobModel_->CHFSubCool
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    )
                );

                const scalarField CHFtotal(CHF*CHFSubCool);

                tDNB =
                    TDNBModel_->TDNB
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    );

                const scalarField MHF
                (
                    MHFModel_->MHF
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    )
                );

                TLeiden =
                    LeidenfrostModel_->TLeid
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    );

                // htc for film boiling
                htcFilmBoiling =
                    filmBoilingModel_->htcFilmBoil
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    );

                // htc for film transition boiling

                // Indicator between CHF (phi = 0) and MHF (phi = 1)
                const scalarField phi
                (
                    min
                    (
                        max
                        (
                            wp_*(Tw - tDNB)/(TLeiden - tDNB),
                            scalar(0)
                        ),
                        scalar(1)
                    )
                );

                Qtb = CHFtotal*(1 - phi) + phi*MHF;

            }


            // Sub-cool boiling Nucleation
            const scalarField N
            (
                nucleationSiteModel_->N
                (
                    liquid,
                    vapor,
                    patchi,
                    Tl,
                    Tsatw,
                    L
                )
            );

            // Bubble departure diameter:
            dDep_ = departureDiamModel_->dDeparture
            (
                liquid,
                vapor,
                patchi,
                Tl,
                Tsatw,
                L
            );

            // Bubble departure frequency:
            const scalarField fDep
            (
                departureFreqModel_->fDeparture
                (
                    liquid,
                    vapor,
                    patchi,
                    dDep_
                )
            );

            // Convective thermal diffusivity for single phase
            alphatConv_ = calcAlphat(alphatConv_);

            // Convective heat transfer area for Sub-cool boiling
            scalarField A1(this->size(), 0);

            const scalarField hewSn(hew.snGrad());

            scalarField alphaFilm(this->size(), 0);

            // Use to identify regimes per face
            labelField regimeTypes(A1.size(), -1);

            forAll(*this, i)
            {
                if (Tw[i] > Tsatw[i])
                {
                    // Sub-cool boiling
                    if (Tw[i] < tDNB[i])
                    {
                        // Sub-cool boiling
                        regimeTypes[i] = regimeType::subcool;

                        // Del Valle & Kenning (1985)
                        const scalar Ja =
                            rhow[i]*Cpw[i]*(Tsatw[i] - Tl[i])
                            /(rhoVaporw[i]*L[i]);

                        const scalar Al =
                            fLiquid[i]*4.8*exp(min(-Ja/80, log(VGREAT)));

                        const scalar A2
                        (
                            min(pi*sqr(dDep_[i])*N[i]*Al/4, scalar(1))
                        );

                        A1[i] = max(1 - A2, scalar(1e-4));

                        // Following Bowring(1962)
                        const scalar A2E
                        (
                            min
                            (
                                pi*sqr(dDep_[i])*N[i]*Al/4,
                                scalar(5)
                            )
                        );

                        // Volumetric mass source in the near wall cell due
                        // to the wall boiling
                        dmdt_[i] =
                            (
                                (1 - relax)*dmdt_[i]
                                + relax*(1.0/6.0)*A2E*dDep_[i]*rhoVaporw[i]
                                * fDep[i]*AbyV_[i]
                            );

                        // Volumetric source in the near wall cell due to
                        // the wall boiling
                        mDotL_[i] = dmdt_[i]*L[i];

                        // Quenching heat transfer coefficient
                        const scalar hQ
                        (
                            2*(alphaw[i]*Cpw[i])*fDep[i]
                            *sqrt
                            (
                                (0.8/max(fDep[i], SMALL))/(pi*alphaw[i]/rhow[i])
                            )
                        );

                        // Quenching heat flux in Sub-cool boiling
                        qq_[i] =
                            (
                                (1 - relax)*qq_[i]
                                + relax*A2*hQ*max(Tw[i] - Tl[i], scalar(0))
                            );

                        this->operator[](i) =
                        (
                            max
                            (
                                A1[i]*alphatConv_[i]
                                + (
                                    (qq_[i] + mDotL_[i]/AbyV_[i])
                                    / max(hewSn[i], scalar(1e-16))
                                )
                                /max(liquidw[i], scalar(1e-8)),
                                1e-8
                            )
                        );
                    }
                    else if (Tw[i] > tDNB[i] && Tw[i] < TLeiden[i])
                    {
                        // transient boiling
                        regimeTypes[i] = regimeType::transient;

                        // No convective heat transfer
                        alphatConv_[i] = 0.0;

                        // transient boiling
                        dmdt_[i] =
                                fLiquid[i]
                                *(
                                    relax*Qtb[i]*AbyV_[i]/L[i]
                                    + (1 - relax)*dmdt_[i]
                                );

                        mDotL_[i] = dmdt_[i]*L[i];


                        // No quenching flux
                        qq_[i] = 0.0;

                        this->operator[](i) =
                        max
                        (
                            (
                                mDotL_[i]/AbyV_[i]
                                /max(hewSn[i], scalar(1e-16))
                            )/max(liquidw[i], scalar(1e-8)),
                            1e-8
                        );

                    }
                    else if (Tw[i] > TLeiden[i])
                    {
                        regimeTypes[i] = regimeType::film; // film boiling

                        // No convective heat transfer
                        alphatConv_[i] = 0.0;

                        // Film boiling
                        dmdt_[i] =
                            fLiquid[i]
                            *(
                                relax*htcFilmBoiling[i]
                                *max(Tw[i] - Tsatw[i], 0)
                                *AbyV_[i]/L[i]
                                + (1 - relax)*dmdt_[i]
                            );


                        mDotL_[i] = dmdt_[i]*L[i];

                        // No quenching flux
                        qq_[i] = 0.0;

                        alphaFilm[i] =
                        (
                            mDotL_[i]/AbyV_[i]/max(hewSn[i], scalar(1e-16))
                        );

                        // alphat is added alphal and multiplied by phase
                        // alphaFilm in the coupled BC. We subtract
                        // alpha and divide by phase to get a net alphaFilm
                        this->operator[](i) =
                            (
                                alphaFilm[i]/max(liquidw[i], scalar(1e-8))
                              - alphaw[i]
                            );
                    }
                }
                else
                {
                    // Tw below Tsat. No boiling single phase convection
                    // Single phase
                    regimeTypes[i] = regimeType::nonBoiling;

                    qq_[i] = 0.0;
                    mDotL_[i] = 0.0;
                    dmdt_[i] = 0.0;

                    // Turbulente thermal diffusivity for single phase.
                    this->operator[](i) =
                    (
                        max
                        (
                            fLiquid[i]*(alphatConv_[i])
                            /max(liquidw[i], scalar(1e-8)),
                            1e-8
                        )
                    );
                }
            }

            if (debug)
            {
                const scalarField qEff
                (
                    liquidw*(*this + alphaw)*hew.snGrad()
                );

                Info<< "alphat for liquid:  " <<  nl << endl;

                Info<< "  alphatl: " << gMin((*this)) << " - "
                    << gMax((*this)) << endl;

                Info<< "  dmdt: " << gMin((dmdt_)) << " - "
                    << gMax((dmdt_)) << endl;

                Info<< "  alphatlEff: " << gMin(liquidw*(*this + alphaw))
                    << " - " << gMax(liquidw*(*this + alphaw)) << endl;

                scalar Qeff = gSum(qEff*patch().magSf());
                Info<< " Effective heat transfer rate to liquid: " << Qeff
                    << endl << nl;

                if (debug & 2)
                {
                    scalarField nSubCools(this->size(), 0);
                    scalarField nTransients(this->size(), 0);
                    scalarField nFilms(this->size(), 0);
                    scalarField nNonBoilings(this->size(), 0);

                    forAll(*this, i)
                    {
                        //faceRegimes[i] = regimeTypes[i];
                        switch (regimeTypes[i])
                        {
                            case regimeType::subcool:
                                nSubCools[i] = 1;
                            break;

                            case regimeType::transient:
                                nTransients[i] = 1;
                            break;

                            case regimeType::film:
                                nFilms[i] = 1;
                            break;

                            case regimeType::nonBoiling:
                                nNonBoilings[i] = 1;
                            break;
                        }
                    }

                    scalar nSubCool(gSum(nSubCools));
                    scalar nTransient(gSum(nTransients));
                    scalar nFilm(gSum(nFilms));
                    scalar nNonBoiling(gSum(nNonBoilings));

                    Info<< "Faces regime :  " <<  nl << endl;

                    Info<< "    sub Cool faces : " << nSubCool << endl;
                    Info<< "    transient faces : " << nTransient << endl;
                    Info<< "    film faces : " << nFilm << endl;
                    Info<< "    non-Boiling faces : " << nNonBoiling << endl;
                    Info<< "    total faces : "
                        << nSubCool + nTransient + nFilm  + nNonBoiling
                        << endl << nl;

                    const scalarField qc
                    (
                        nNonBoilings*fLiquid*A1*(alphatConv_ + alphaw)
                        *hew.snGrad()
                    );

                    scalar Qc = gSum(qc*patch().magSf());
                    Info<< " Convective heat transfer: " << Qc << endl;

                    const scalarField qFilm
                    (
                        relax*fLiquid*nFilms*htcFilmBoiling*(Tw - Tsatw)
                    );

                    scalar QFilm = gSum(qFilm*patch().magSf());
                    Info<< " Film boiling heat transfer: " << QFilm << endl;

                    Info<< " Htc Film Boiling coeff: "
                        << gMin(nFilms*htcFilmBoiling)
                        << " - "
                        << gMax(nFilms*htcFilmBoiling) << endl;

                    scalar Qtbtot =
                        gSum(fLiquid*nTransients*Qtb*patch().magSf());
                    Info<< " Transient boiling heat transfer:" << Qtbtot
                        << endl;


                    Info<< " TDNB: " << gMin(tDNB) << " - " << gMax(tDNB)
                        << endl;

                    const scalarField qSubCool
                    (
                        fLiquid*nSubCools*
                        (
                            A1*alphatConv_*hew.snGrad()
                            + qe() + qq()
                        )
                    );


                    scalar QsubCool = gSum(qSubCool*patch().magSf());

                    Info<< " Sub Cool boiling heat transfer: " << QsubCool
                        << endl;

                    Info<< "  N: " << gMin(nSubCools*N) << " - "
                        << gMax(nSubCools*N) << endl;
                    Info<< "  dDep: " << gMin(nSubCools*dDep_) << " - "
                        << gMax(nSubCools*dDep_) << endl;
                    Info<< "  fDep: " << gMin(nSubCools*fDep) << " - "
                        << gMax(nSubCools*fDep) << endl;
                    Info<< "  A1: " << gMin(nSubCools*A1) << " - "
                        << gMax(nSubCools*A1) << endl;

                    Info<< nl;
                }
            }
        }
        break;
        default:
        {
            FatalErrorInFunction
                << "Unknown phase type. Valid types are: "
                << phaseTypeNames_ << nl << exit(FatalError);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallBoilingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("phaseType", phaseTypeNames_[phaseType_]);

    relax_->writeData(os);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            os.beginBlock("partitioningModel");
            partitioningModel_->write(os);
            os.endBlock();

            if (filmBoilingModel_)
            {
                os.beginBlock("filmBoilingModel");
                filmBoilingModel_->write(os);
                os.endBlock();
            }

            if (LeidenfrostModel_)
            {
                os.beginBlock("LeidenfrostModel");
                LeidenfrostModel_->write(os);
                os.endBlock();
            }

            break;
        }
        case liquidPhase:
        {
            os.beginBlock("partitioningModel");
            partitioningModel_->write(os);
            os.endBlock();

            os.beginBlock("nucleationSiteModel");
            nucleationSiteModel_->write(os);
            os.endBlock();

            os.beginBlock("departureDiamModel");
            departureDiamModel_->write(os);
            os.endBlock();

            os.beginBlock("departureFreqModel");
            departureFreqModel_->write(os);
            os.endBlock();

            if (filmBoilingModel_)
            {
                os.beginBlock("filmBoilingModel");
                filmBoilingModel_->write(os);
                os.endBlock();
            }

            if (LeidenfrostModel_)
            {
                os.beginBlock("LeidenfrostModel");
                LeidenfrostModel_->write(os);
                os.endBlock();
            }

            if (CHFModel_)
            {
                os.beginBlock("CHFModel");
                CHFModel_->write(os);
                os.endBlock();
            }

            if (CHFSoobModel_)
            {
                os.beginBlock("CHFSubCoolModel");
                CHFSoobModel_->write(os);
                os.endBlock();
            }

            if (MHFModel_)
            {
                os.beginBlock("MHFModel");
                MHFModel_->write(os);
                os.endBlock();
            }

            if (TDNBModel_)
            {
                os.beginBlock("TDNBModel");
                TDNBModel_->write(os);
                os.endBlock();
            }

            os.writeEntry("K", K_);
            os.writeEntry("wp", wp_);
            break;
        }
    }

    os.writeEntry("otherPhase", otherPhaseName_);

    dmdt_.writeEntry("dmdt", os);
    dDep_.writeEntry("dDep", os);
    qq_.writeEntry("qQuenching", os);
    alphatConv_.writeEntry("alphatConv", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatWallBoilingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
