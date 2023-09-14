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

\*---------------------------------------------------------------------------*/

#include "turbulentTemperatureRadCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "basicThermo.H"
#include "IOField.H"
#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField&
turbulentTemperatureRadCoupledMixedFvPatchScalarField::getOrCreateField
(
    const word& fieldName
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    auto* ptr = mesh.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
        mesh.objectRegistry::store(ptr);
    }

    return *ptr;
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::storeHTCFields
(
    const word& prefix,
    const scalarField& shtc,
    const scalarField& shtcPatch
)
const
{
    volScalarField& htc =
        getOrCreateField(IOobject::scopedName(prefix, "htc"));
    htc.boundaryFieldRef()[patch().index()] = shtc;

    volScalarField& htcPatch =
        getOrCreateField(IOobject::scopedName(prefix, "htcPatch"));
    htcPatch.boundaryFieldRef()[patch().index()] = shtcPatch;
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::writeFileHeader
(
    Ostream& os
)
{
    writeCommented(os, "Time");
    writeTabbed(os, "Q_[W]");
    writeTabbed(os, "q_[W/m^2]");
    writeTabbed(os, "HTCavg_[W/m^2/K]");
    writeTabbed(os, "patchHTCavg_[W/m^2/K]");
    writeTabbed(os, "TpMin_[K]");
    writeTabbed(os, "TpMax_[K]");
    writeTabbed(os, "TpAvg_[K]");
    writeTabbed(os, "TpNbrMin_[K]");
    writeTabbed(os, "TpNbrMax_[K]");
    writeTabbed(os, "TpNbrAvg_[K]");

    os  << endl;

    writtenHeader_ = true;
    updateHeader_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),  // default method (fluidThermo)
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this
    ),
    functionObjects::writeFile
    (
        db(),
        "turbulentTemperatureRadCoupledMixed",
        "undefined",
        false
    ),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr"),
    logInterval_(-1),
    executionIndex_(0),
    thermalInertia_(false),
    verbose_(false),
    prefix_(word::null)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
    this->source() = 0.0;
}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        psf
    ),
    functionObjects::writeFile(psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    thicknessLayer_(psf.thicknessLayer_.clone(p.patch())),
    kappaLayers_(psf.kappaLayers_),
    kappaLayer_(psf.kappaLayer_.clone(p.patch())),
    logInterval_(psf.logInterval_),
    executionIndex_(psf.executionIndex_),
    thermalInertia_(psf.thermalInertia_),
    verbose_(psf.verbose_),
    prefix_(psf.prefix_)
{}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        dict
    ),
    functionObjects::writeFile
    (
        db(),
        "turbulentTemperatureRadCoupledMixed",
        patch().name(),
        false
    ),
    TnbrName_(dict.getOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.getOrDefault<word>("qrNbr", "none")),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    logInterval_(dict.getOrDefault<scalar>("logInterval", -1)),
    executionIndex_(0),
    thermalInertia_(dict.getOrDefault<Switch>("thermalInertia", false)),
    verbose_(dict.getOrDefault<bool>("verbose", false)),
    prefix_(dict.getOrDefault<word>("prefix", "multiWorld"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    //const auto* eptr = dict.findEntry("thicknessLayers", keyType::LITERAL);
    //if (eptr)
    //{
    //    // Detect either a list (parsed as a scalarList) or
    //    // a single entry (parsed as a PatchFunction1) or
    //
    //    if
    //    (
    //        eptr->isStream()
    //     && eptr->stream().peek().isPunctuation(token::BEGIN_LIST)
    //    )
    //    {
    //        // Backwards compatibility
    //        thicknessLayers_ = dict.get<scalarList>("thicknessLayers");
    //        kappaLayers_ = dict.get<scalarList>("kappaLayers");
    //
    //        if (thicknessLayers_.size() != kappaLayers_.size())
    //        {
    //            FatalIOErrorInFunction(dict) << "Inconstent sizes :"
    //                << "thicknessLayers:" << thicknessLayers_
    //                << "kappaLayers:" << kappaLayers_
    //                << exit(FatalIOError);
    //        }
    //    }
    //    else
    //    {
    //        thicknessLayer_ = PatchFunction1<scalar>::New
    //        (
    //            p.patch(),
    //            "thicknessLayers",
    //            dict
    //        );
    //        kappaLayer_ = PatchFunction1<scalar>::New
    //        (
    //            p.patch(),
    //            "kappaLayers",
    //            dict
    //        );
    //    }
    //}

    // Read list of layers
    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }
    // Read single additional layer as PatchFunction1
    thicknessLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "thicknessLayer",
        dict
    );
    if (thicknessLayer_)
    {
        kappaLayer_ = PatchFunction1<scalar>::New
        (
            p.patch(),
            "kappaLayer",
            dict
        );
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
        refGrad() = Zero;
        valueFraction() = 1.0;
    }

    bool boolVal(false);
    if (dict.readIfPresent("useImplicit", boolVal))
    {
        this->useImplicit(boolVal);
    }
    if (dict.found("source"))
    {
        source() = scalarField("source", dict, p.size());
    }
    else
    {
        source() = 0.0;
    }

    writeFile::read(dict);
}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        psf
    ),
    functionObjects::writeFile(psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    thicknessLayer_(psf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(psf.kappaLayers_),
    kappaLayer_(psf.kappaLayer_.clone(patch().patch())),
    logInterval_(psf.logInterval_),
    executionIndex_(psf.executionIndex_),
    thermalInertia_(psf.thermalInertia_),
    verbose_(psf.verbose_),
    prefix_(psf.prefix_)
{}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf
)
:
    mixedFvPatchScalarField(psf),
    temperatureCoupledBase(patch(), psf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), psf.internalField()),
        *this,
        psf
    ),
    functionObjects::writeFile(psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    thicknessLayer_(psf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(psf.kappaLayers_),
    kappaLayer_(psf.kappaLayer_.clone(patch().patch())),
    logInterval_(psf.logInterval_),
    executionIndex_(psf.executionIndex_),
    thermalInertia_(psf.thermalInertia_),
    verbose_(psf.verbose_),
    prefix_(psf.prefix_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureRadCoupledMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
    if (thicknessLayer_)
    {
        thicknessLayer_().autoMap(mapper);
        kappaLayer_().autoMap(mapper);
    }
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& tiptf =
        refCast
        <
            const turbulentTemperatureRadCoupledMixedFvPatchScalarField
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
    if (thicknessLayer_)
    {
        thicknessLayer_().rmap(tiptf.thicknessLayer_(), addr);
        kappaLayer_().rmap(tiptf.kappaLayer_(), addr);
    }
}


tmp<Foam::scalarField>
turbulentTemperatureRadCoupledMixedFvPatchScalarField::kappa
(
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));

    return tk;
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const label patchi = patch().index();
    const mappedPatchBase& mpp =
        mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );

    const scalarField Tc(patchInternalField());
    const scalarField& Tp = *this;

    const scalarField kappaTp(kappa(Tp));
    const scalarField KDelta(kappaTp*patch().deltaCoeffs());


    scalarField TcNbr;
    scalarField TpNbr;
    scalarField KDeltaNbr;

    if (mpp.sameWorld())
    {
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const label samplePatchi = mpp.samplePolyPatch().index();
        const fvPatch& nbrPatch =
            refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

        const auto& nbrField = refCast
                <const turbulentTemperatureRadCoupledMixedFvPatchScalarField>
                (
                    nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
                );

        // Swap to obtain full local values of neighbour K*delta
        TcNbr = nbrField.patchInternalField();
        TpNbr = nbrField;
        KDeltaNbr = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering.
        TcNbr = patchInternalField();
        TpNbr = Tp;
        KDeltaNbr = KDelta;
    }
    distribute(this->internalField().name() + "_value", TcNbr);
    distribute(this->internalField().name() + "_patchValue", TpNbr);
    distribute(this->internalField().name() + "_weights", KDeltaNbr);

    scalarField KDeltaC(this->size(), GREAT);
    if (thicknessLayer_ || thicknessLayers_.size())
    {
        // Harmonic averaging
        {
            KDeltaC = 0.0;

            if (thicknessLayer_)
            {
                const scalar t = db().time().timeOutputValue();
                KDeltaC +=
                     thicknessLayer_().value(t)
                    /kappaLayer_().value(t);

            }
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    KDeltaC += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
                }
            }
            KDeltaC = 1.0/(KDeltaC + SMALL);
        }
    }

    scalarField alpha(kappaTp*(1 + KDeltaNbr/KDeltaC)*patch().deltaCoeffs());

    scalarField qr(Tp.size(), Zero);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField>(qrName_);
    }

    scalarField qrNbr(Tp.size(), Zero);
    if (qrNbrName_ != "none")
    {
        if (mpp.sameWorld())
        {
            const polyMesh& nbrMesh = mpp.sampleMesh();
            const label samplePatchi = mpp.samplePolyPatch().index();
            const fvPatch& nbrPatch =
                refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];
            qrNbr = nbrPatch.lookupPatchField<volScalarField>(qrNbrName_);
        }
        else
        {
            qrNbr = patch().lookupPatchField<volScalarField>(qrNbrName_);
        }
        distribute(qrNbrName_, qrNbr);
    }

    // inertia therm
    if (thermalInertia_ && !mpp.sameWorld())
    {
        FatalErrorInFunction
            << "thermalInertia not supported in combination with multi-world"
            << exit(FatalError);
    }
    if (thermalInertia_)
    {
        const scalar dt = mesh.time().deltaTValue();
        scalarField mCpDtNbr;

        {
            const polyMesh& nbrMesh = mpp.sampleMesh();

            const basicThermo* thermo =
                nbrMesh.findObject<basicThermo>(basicThermo::dictName);

            if (thermo)
            {
                const label samplePatchi = mpp.samplePolyPatch().index();
                const fvPatch& nbrPatch =
                    refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];
                const scalarField& ppn =
                    thermo->p().boundaryField()[samplePatchi];
                const scalarField& Tpn =
                    thermo->T().boundaryField()[samplePatchi];

                mCpDtNbr =
                (
                    thermo->Cp(ppn, Tpn, samplePatchi)
                  * thermo->rho(samplePatchi)
                  / nbrPatch.deltaCoeffs()/dt
                );

                mpp.distribute(mCpDtNbr);
            }
            else
            {
                mCpDtNbr.setSize(Tp.size(), Zero);
            }
        }

        scalarField mCpDt;

        // Local inertia therm
        {
            const basicThermo* thermo =
                mesh.findObject<basicThermo>(basicThermo::dictName);

            if (thermo)
            {
                const scalarField& pp = thermo->p().boundaryField()[patchi];

                mCpDt =
                (
                    thermo->Cp(pp, Tp, patchi)
                  * thermo->rho(patchi)
                  / patch().deltaCoeffs()/dt
                );
            }
            else
            {
                // Issue warning?
                mCpDt.setSize(Tp.size(), Zero);
            }
        }

        const volScalarField& T =
            this->db().lookupObject<volScalarField>
            (
                this->internalField().name()
            );

        const fvPatchField<scalar>& TpOld = T.oldTime().boundaryField()[patchi];

        scalarField alpha(KDeltaNbr + mCpDt + mCpDtNbr);

        valueFraction() = alpha/(alpha + KDelta);
        scalarField c(KDeltaNbr*TcNbr + (mCpDt + mCpDtNbr)*TpOld);
        refValue() = c/alpha;
        refGrad() = (qr + qrNbr)/kappaTp;
    }
    else
    {

        valueFraction() = KDeltaNbr/(KDeltaNbr + alpha);
        refValue() = TcNbr;
        refGrad() = (qr + qrNbr)/kappaTp;
    }

    source() = Zero;

    // If useImplicit is true we need the source term associated with this BC
    if (this->useImplicit())
    {
        source() =
            alphaSfDelta()*
            (
                valueFraction()*deltaH()
              + (qr + qrNbr)/beta()
            );
    }

    mixedFvPatchScalarField::updateCoeffs();


    if (verbose_)
    {
        // Calculate heat-transfer rate and heat flux
        const scalar Q = gSum(kappaTp*patch().magSf()*snGrad());
        const scalar magSf = gSum(patch().magSf());
        const scalar q = Q/max(magSf, SMALL);


        // Calculate heat-transfer coeff based on the first definition
        // [W/m^2] = [W/m/K K * 1/m]
        const scalarField qField
        (
            kappaTp*snGrad()
        );
        const scalarField deltaT(TcNbr - Tc);
        scalarField htc(deltaT.size(), Zero);

        forAll(deltaT, i)
        {
            if (mag(deltaT[i]) > SMALL)
            {
                htc[i] = qField[i]/deltaT[i];
            }
        }
        const scalar aveHtc = gSum(htc*patch().magSf())/max(magSf, SMALL);


        // Calculate heat-transfer coeff based on the second definition
        const scalarField deltaTPatch(TpNbr - Tp);
        scalarField htcPatch(deltaTPatch.size(), Zero);

        forAll(deltaTPatch, i)
        {
            if (mag(deltaTPatch[i]) > SMALL)
            {
                htcPatch[i] = qField[i]/deltaTPatch[i];
            }
        }
        const scalar aveHtcPatch =
            gSum(htcPatch*patch().magSf())/max(magSf, SMALL);

        // Calculate various averages of temperature
        const scalarMinMax TpMinMax = gMinMax(Tp);
        const scalar TpAvg = gAverage(Tp);
        const scalarMinMax TpNbrMinMax = gMinMax(TpNbr);
        const scalar TpNbrAvg = gAverage(TpNbr);


        Info<< nl
            << patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :" << nl
            << " Heat transfer rate [W]:" << Q << nl
            << " Area [m^2]:" << magSf << nl
            << " Heat flux [W/m^2]:" << q << nl
            << " Area-averaged heat-transfer coefficient [W/m^2/K]:"
            << aveHtc << nl
            << " Area-averaged patch heat-transfer coefficient [W/m^2/K]:"
            << aveHtcPatch << nl
            << " Wall temperature [K]"
            << " min:" << TpMinMax.min()
            << " max:" << TpMinMax.max()
            << " avg:" << TpAvg << nl
            << " Neighbour wall temperature [K]"
            << " min:" << TpNbrMinMax.min()
            << " max:" << TpNbrMinMax.max()
            << " avg:" << TpNbrAvg
            << nl << endl;


        // Handle data for file output
        if (canResetFile())
        {
            resetFile(patch().name());
        }

        if (canWriteHeader())
        {
            writeFileHeader(file());
        }

        if (canWriteToFile() && writeFile())
        {
            file()
                << db().time().timeOutputValue() << token::TAB
                << Q << token::TAB
                << q << token::TAB
                << aveHtc << token::TAB
                << aveHtcPatch << token::TAB
                << TpMinMax.min() << token::TAB
                << TpMinMax.max() << token::TAB
                << TpAvg << token::TAB
                << TpNbrMinMax.min() << token::TAB
                << TpNbrMinMax.max() << token::TAB
                << TpNbrAvg << token::TAB
                << endl;
        }


        // Store htc fields as patch fields of a volScalarField
        storeHTCFields(prefix_, htc, htcPatch);
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& m,
    const label iMatrix,
    const direction cmpt
)
{
    FatalErrorInFunction
        << "This T BC does not support energy coupling "
        << "It is implemented on he field "
        << abort(FatalError);
}


tmp<Field<scalar>> turbulentTemperatureRadCoupledMixedFvPatchScalarField::coeffs
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    FatalErrorInFunction
        << "This BC does not support energy coupling "
        << "Use compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << abort(FatalError);

    return tmp<Field<scalar>>(new Field<scalar>());
}


tmp<scalarField>
turbulentTemperatureRadCoupledMixedFvPatchScalarField::alphaSfDelta() const
{
    return (alpha(*this)*patch().deltaCoeffs()*patch().magSf());
}


tmp<scalarField> turbulentTemperatureRadCoupledMixedFvPatchScalarField::
beta() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    if (!mpp.sameWorld())
    {
        FatalErrorInFunction
            << "coupled energy not supported in combination with multi-world"
            << exit(FatalError);
    }

    const label samplePatchi = mpp.samplePolyPatch().index();
    const polyMesh& nbrMesh = mpp.sampleMesh();

    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    const turbulentTemperatureRadCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const turbulentTemperatureRadCoupledMixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
            );

    // Swap to obtain full local values of neighbour internal field
    scalarField TcNbr(nbrField.patchInternalField());
    mpp.distribute(TcNbr);

    scalarField alphaDeltaNbr
    (
        nbrField.alpha(TcNbr)*nbrPatch.deltaCoeffs()
    );
    mpp.distribute(alphaDeltaNbr);

    scalarField alphaDelta
    (
         alpha(*this)*patch().deltaCoeffs()
    );

    return (alphaDeltaNbr + alphaDelta);
}


tmp<scalarField> turbulentTemperatureRadCoupledMixedFvPatchScalarField::
deltaH() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    if (!mpp.sameWorld())
    {
        FatalErrorInFunction
            << "coupled energy not supported in combination with multi-world"
            << exit(FatalError);
    }

    const polyMesh& nbrMesh = mpp.sampleMesh();

    const basicThermo* nbrThermo =
        nbrMesh.cfindObject<basicThermo>(basicThermo::dictName);

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    const basicThermo* localThermo =
        mesh.cfindObject<basicThermo>(basicThermo::dictName);


    if (nbrThermo && localThermo)
    {
        const label patchi = patch().index();
        const scalarField& pp = localThermo->p().boundaryField()[patchi];
        const scalarField& Tp = *this;

        const mappedPatchBase& mpp =
            refCast<const mappedPatchBase>(patch().patch());

        const label patchiNrb = mpp.samplePolyPatch().index();

        const scalarField& ppNbr = nbrThermo->p().boundaryField()[patchiNrb];
        //const scalarField& TpNbr = nbrThermo->T().boundaryField()[patchiNrb];

        // Use this Tp to evaluate he jump as this is updated while doing
        // updateCoeffs on boundary fields which set T values on patches
        // then non consistent Tp and Tpnbr could be used from different
        // updated values (specially when T changes drastically between time
        // steps/
        return
        (
          -  localThermo->he(pp, Tp, patchi)
          +  nbrThermo->he(ppNbr, Tp, patchiNrb)
        );
    }
    else
    {
        FatalErrorInFunction
            << "Can't find thermos on mapped patch "
            << " method, but thermo package not available"
            << exit(FatalError);
    }

    return tmp<scalarField>::New(patch().size(), Zero);
}


bool turbulentTemperatureRadCoupledMixedFvPatchScalarField::writeFile()
{
    if (!verbose_ || (logInterval_ <= 0))
    {
        return false;
    }

    const auto& time = patch().boundaryMesh().mesh().time();

    const scalar t = time.timeOutputValue();
    const scalar ts = time.startTime().value();
    const scalar deltaT = time.deltaTValue();

    const label executionIndex = label
    (
        (
            (t - ts)
          + 0.5*deltaT
        )
        /logInterval_
    );

    bool write = false;
    if (executionIndex > executionIndex_)
    {
        executionIndex_ = executionIndex;
        write = true;
    }

    return write;
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);

    os.writeEntryIfDifferent<word>("Tnbr", "T", TnbrName_);
    os.writeEntryIfDifferent<word>("qrNbr", "none", qrNbrName_);
    os.writeEntryIfDifferent<word>("qr", "none", qrName_);
    os.writeEntry<scalar>("logInterval", logInterval_);

    if (thermalInertia_)
    {
        os.writeEntry("thermalInertia", thermalInertia_);
    }
    os.writeEntryIfDifferent<bool>("verbose", false, verbose_);
    os.writeEntryIfDifferent<word>("prefix", "multiWorld", prefix_);

    if (thicknessLayer_)
    {
        thicknessLayer_().writeData(os);
        kappaLayer_().writeData(os);
    }
    if (thicknessLayers_.size())
    {
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }

    // Write writeFile entries
    os.writeEntry<label>("writePrecision", writePrecision_);
    os.writeEntry<bool>("updateHeader", updateHeader_);
    os.writeEntry<bool>("writeToFile", writeToFile_);
    os.writeEntry<bool>("useUserTime", useUserTime_);

    temperatureCoupledBase::write(os);
    mappedPatchFieldBase<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureRadCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
