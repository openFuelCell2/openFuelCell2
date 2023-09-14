/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentFluidThermoModel.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    qs_(p.size()),
    solidDict_(),
    solidPtr_(nullptr),
    qrPrevious_(p.size()),
    qrRelaxation_(1),
    qrName_("undefined-qr")
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase(p.patch(), ptf),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_, mapper),
    qs_(ptf.qs_, mapper),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    qrPrevious_(ptf.qrPrevious_, mapper),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), NEARESTPATCHFACE, dict),
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(dict.getOrDefault("baffleActivated", true)),
    thickness_(),
    qs_(p.size(), Zero),
    solidDict_(dict),
    solidPtr_(),
    qrPrevious_(p.size(), Zero),
    qrRelaxation_
    (
        dict.getOrDefaultCompat("qrRelaxation", {{"relaxation", 1712}}, 1)
    ),
    qrName_(dict.getOrDefault<word>("qr", "none"))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dict, p.size());
    }

    if (dict.found("qs"))
    {
        qs_ = scalarField("qs", dict, p.size());
    }

    if (dict.found("qrPrevious"))
    {
        qrPrevious_ = scalarField("qrPrevious", dict, p.size());
    }

    if (baffleActivated_ && this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }

}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    qrPrevious_(ptf.qrPrevious_),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    qrPrevious_(ptf.qrPrevious_),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class solidType>
bool thermalBaffle1DFvPatchScalarField<solidType>::owner() const
{
    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    return (patchi < nbrPatchi);
}


template<class solidType>
const solidType& thermalBaffle1DFvPatchScalarField<solidType>::solid() const
{
    if (this->owner())
    {
        if (!solidPtr_)
        {
            solidPtr_.reset(new solidType(solidDict_));
        }
        return *solidPtr_;
    }
    else
    {
        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField>(TName_)
        );

        return nbrField.solid();
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::
baffleThickness() const
{
    if (this->owner())
    {
        if (thickness_.size() != patch().size())
        {
            FatalIOErrorInFunction(solidDict_)
                << "Field thickness has not been specified"
                   " for patch " << this->patch().name()
                << exit(FatalIOError);
        }

        return thickness_;
    }
    else
    {
        const mapDistribute& mapDist = this->mappedPatchBase::map();

        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];
        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField>(TName_)
        );

        tmp<scalarField> tthickness
        (
            new scalarField(nbrField.baffleThickness())
        );
        scalarField& thickness = tthickness.ref();
        mapDist.distribute(thickness);
        return tthickness;
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::qs() const
{
    if (this->owner())
    {
         return qs_;
    }
    else
    {
        const mapDistribute& mapDist = this->mappedPatchBase::map();

        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField>(TName_)
        );

        tmp<scalarField> tqs(new scalarField(nbrField.qs()));
        scalarField& qs = tqs.ref();
        mapDist.distribute(qs);
        return tqs;
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mappedPatchBase::clearOut();

    mixedFvPatchScalarField::autoMap(m);

    if (this->owner())
    {
        thickness_.autoMap(m);
        qs_.autoMap(m);
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.rmap(tiptf.thickness_, addr);
        qs_.rmap(tiptf.qs_, addr);
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const mapDistribute& mapDist = this->mappedPatchBase::map();

    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    if (baffleActivated_)
    {
        const fvPatch& nbrPatch = patch().boundaryMesh()[nbrPatchi];

        const compressible::turbulenceModel& turbModel =
            db().template lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        // local properties
        const scalarField kappaw(turbModel.kappaEff(patchi));

        const fvPatchScalarField& Tp =
            patch().template lookupPatchField<volScalarField>(TName_);


        scalarField qr(Tp.size(), Zero);

        if (qrName_ != "none")
        {
            qr = patch().template lookupPatchField<volScalarField>(qrName_);

            qr = qrRelaxation_*qr + (1.0 - qrRelaxation_)*qrPrevious_;
            qrPrevious_ = qr;
        }

        tmp<scalarField> Ti = patchInternalField();

        scalarField myKDelta(patch().deltaCoeffs()*kappaw);

        // nrb properties
        scalarField nbrTp =
            turbModel.transport().T().boundaryField()[nbrPatchi];
        mapDist.distribute(nbrTp);

        // solid properties
        scalarField kappas(patch().size(), Zero);
        forAll(kappas, i)
        {
            kappas[i] = solid().kappa(0.0, (Tp[i] + nbrTp[i])/2.0);
        }

        scalarField KDeltaSolid(kappas/baffleThickness());

        scalarField alpha(KDeltaSolid - qr/Tp);

        valueFraction() = alpha/(alpha + myKDelta);

        refValue() = (KDeltaSolid*nbrTp + qs()/2.0)/alpha;

        if (debug)
        {
            scalar Q = gAverage(kappaw*snGrad());
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrPatch.name() << ':'
                << this->internalField().name() << " :"
                << " heat[W]:" << Q
                << " walltemperature "
                << " min:" << gMin(*this)
                << " max:" << gMax(*this)
                << " avg:" << gAverage(*this)
                << endl;
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}

template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);
    mappedPatchBase::write(os);

    if (this->owner())
    {
        baffleThickness()().writeEntry("thickness", os);
        qs()().writeEntry("qs", os);
        solid().write(os);
    }

    qrPrevious_.writeEntry("qrPrevious", os);
    os.writeEntry("qr", qrName_);
    os.writeEntry("qrRelaxation", qrRelaxation_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
