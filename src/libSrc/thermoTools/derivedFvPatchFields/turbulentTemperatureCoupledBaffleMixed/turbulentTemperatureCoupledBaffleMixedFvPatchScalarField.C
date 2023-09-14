/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
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
    TnbrName_("undefined-Tnbr")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        ptf
    ),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    thicknessLayer_(ptf.thicknessLayer_.clone(p.patch())),
    kappaLayers_(ptf.kappaLayers_),
    kappaLayer_(ptf.kappaLayer_.clone(p.patch()))
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
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
    TnbrName_(dict.get<word>("Tnbr"))
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

    WarningInFunction
        << "This BC has been superseded by "
        << "compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << endl;

    // Read list of layers
    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }
    // Read single additional PatchFunction1
    thicknessLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "thicknessLayer",
        dict
    );
    kappaLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "kappaLayer",
        dict
    );


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

// This blocks (crashes) with more than two worlds!
//
///    // Store patch value as initial guess when running in database mode
///    mappedPatchFieldBase<scalar>::initRetrieveField
///    (
///        this->internalField().name(),
///        *this
///    );
///    mappedPatchFieldBase<scalar>::initRetrieveField
///    (
///        this->internalField().name() + "_weights",
///        this->patch().deltaCoeffs()
///    );
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch()))
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& wtcsf
)
:
    mixedFvPatchScalarField(wtcsf),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), wtcsf.internalField()),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::autoMap
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


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& tiptf =
        refCast
        <
            const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
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
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::kappa
(
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));

    // Optionally modify with explicit resistance
    if (thicknessLayer_ || thicknessLayers_.size())
    {
        scalarField KDelta(tk*patch().deltaCoeffs());

        // Harmonic averaging of kappa*deltaCoeffs
        {
            KDelta = 1.0/KDelta;
            if (thicknessLayer_)
            {
                const scalar t = db().time().timeOutputValue();
                KDelta +=
                    thicknessLayer_().value(t)
                   /kappaLayer_().value(t);
            }
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    KDelta += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
                }
            }
            KDelta = 1.0/KDelta;
        }

        // Update kappa from KDelta
        tk = KDelta/patch().deltaCoeffs();
    }

    return tk;
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );

    const scalarField& Tp = *this;
    const scalarField kappaTp(kappa(Tp));
    const tmp<scalarField> myKDelta = kappaTp*patch().deltaCoeffs();


    scalarField nbrIntFld;
    scalarField nbrKDelta;
    if (mpp.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
        const label nbrPatchID = mpp.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField&
        nbrField =
        refCast
        <
        const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
        );

        // Swap to obtain full local values of neighbour K*delta
        nbrIntFld = nbrField.patchInternalField();
        nbrKDelta = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering.
        nbrIntFld = patchInternalField();
        nbrKDelta = myKDelta.ref();
    }
    distribute(this->internalField().name() + "_value", nbrIntFld);
    distribute(this->internalField().name() + "_weights", nbrKDelta);


    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    this->refValue() = nbrIntFld;
    this->refGrad() = Zero;
    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta());

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappaTp*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& m,
    const label iMatrix,
    const direction cmpt
)
{
    FatalErrorInFunction
        << "This BC does not support energy coupling "
        << "Use compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << abort(FatalError);
}


tmp<Field<scalar>>
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::coeffs
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
    /*
    const label index(this->patch().index());

    const label nSubFaces(matrix.lduMesh().cellBoundMap()[mat][index].size());

    Field<scalar> mapCoeffs(nSubFaces, Zero);

    label subFaceI = 0;
    forAll(*this, faceI)
    {

    }
    */
    return tmp<Field<scalar>>(new Field<scalar>());
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntry("Tnbr", TnbrName_);
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
    temperatureCoupledBase::write(os);
    mappedPatchFieldBase<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
