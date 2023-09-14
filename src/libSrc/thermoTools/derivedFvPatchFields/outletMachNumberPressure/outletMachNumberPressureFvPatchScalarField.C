/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "outletMachNumberPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletMachNumberPressureFvPatchScalarField::
outletMachNumberPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    M_(1),
    pBack_(0.0),
    c1_(0.0),
    A1_(0.0),
    phiName_("phi"),
    rhoName_("rho"),
    UName_("U"),
    choked_(false),
    relax_(0.0)
{}


Foam::outletMachNumberPressureFvPatchScalarField::
outletMachNumberPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    M_(dict.getOrDefault<scalar>("M", 0)),
    pBack_(dict.get<scalar>("pBack")),
    c1_(dict.getOrDefault<scalar>("c1", 0)),
    A1_(dict.getOrDefault<scalar>("A1", 0)),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    UName_(dict.getOrDefault<word>("U", "U")),
    choked_(dict.get<Switch>("choked")),
    relax_(dict.getOrDefault<scalar>("relax", 0))
{}


Foam::outletMachNumberPressureFvPatchScalarField::
outletMachNumberPressureFvPatchScalarField
(
    const outletMachNumberPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    M_(ptf.M_),
    pBack_(ptf.pBack_),
    c1_(ptf.c1_),
    A1_(ptf.A1_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    UName_(ptf.UName_),
    choked_(ptf.choked_),
    relax_(ptf.relax_)
{}


Foam::outletMachNumberPressureFvPatchScalarField::
outletMachNumberPressureFvPatchScalarField
(
    const outletMachNumberPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    M_(tppsf.M_),
    pBack_(tppsf.pBack_),
    c1_(tppsf.c1_),
    A1_(tppsf.A1_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    UName_(tppsf.UName_),
    choked_(tppsf.choked_),
    relax_(tppsf.relax_)
{}


Foam::outletMachNumberPressureFvPatchScalarField::
outletMachNumberPressureFvPatchScalarField
(
    const outletMachNumberPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    M_(tppsf.M_),
    pBack_(tppsf.pBack_),
    c1_(tppsf.c1_),
    A1_(tppsf.A1_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    UName_(tppsf.UName_),
    choked_(tppsf.choked_),
    relax_(tppsf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletMachNumberPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

     const volScalarField& p =
        db().lookupObject<volScalarField>
        (
            this->internalField().name()
        );

    const label patchi = patch().index();

    const scalarField pb(p.oldTime().boundaryField()[patchi]);

    const auto& phi = patch().lookupPatchField<surfaceScalarField>(phiName_);

    // Calculate the current mass flow rate
    if (phi.internalField().dimensions() != dimMass/dimTime)
    {
        FatalErrorInFunction
            <<"phi is not a mass flux." << exit(FatalError);
    }

    const fluidThermo* thermoPtr =
        db().findObject<fluidThermo>(basicThermo::dictName);

    const volVectorField& U = db().lookupObject<volVectorField>(UName_);

    vectorField Ub(U.boundaryField()[patchi]);
    const vectorField UbOld(U.oldTime().boundaryField()[patchi]);

    // relax U
    Ub = relax_*UbOld + (1 - relax_)*Ub;

    const scalarField gamma(thermoPtr->gamma()().boundaryField()[patchi]);

    const auto& rho = patch().lookupPatchField<volScalarField>(rhoName_);

    const scalarField Mb(mag(Ub)/sqrt(gamma*pb/rho));

    const scalarField ptot
    (
        pb*(pow(1 + (gamma-1)/2*sqr(gAverage(Mb)), gamma/(gamma-1)))
    );

    scalarField M(patch().size(), 1.0);

    if (choked_)
    {
        if (M_ > 0.0)
        {
            M = M_;
        }
        else
        {
            FatalErrorInFunction <<" Mach number is lower than zero" << endl
                << "Pelase specify M in the dictionary"
                << exit(FatalError);
        }
    }
    else
    {

        if (A1_ == 0.0 || c1_ == 0.0)
        {
             FatalErrorInFunction <<" Please enter values for A1 and c1" << endl
                << exit(FatalError);
        }

        const scalarField r(pBack_/ptot);
        const scalar area = gSum(mag(patch().Sf()));
        M =
            A1_/(c1_*area)
           *sqrt(2/(gamma-1)*(pow(r, 2/gamma) - pow(r, (gamma+1)/gamma)));

        forAll(M, i)
        {
            if (M[i] < 0 || r[i] >=1)
            {
                WarningInFunction <<" Mach number is lower than zero" << endl
                    << "or pBack/ptot ratio is larger then one"
                    << endl;
            }
        }
    }

    scalarField pbNew
    (
        ptot/(pow(1 + (gamma-1)/2*sqr(gAverage(M)), gamma/(gamma-1)))
    );

    // relax pressure
    pbNew = relax_*pb + (1 -relax_)*pbNew;

    operator==(pbNew);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::outletMachNumberPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("pBack", pBack_);
    os.writeEntryIfDifferent<scalar>("c1", 0, c1_);
    os.writeEntryIfDifferent<scalar>("A1", 0, A1_);
    os.writeEntry("choked", choked_);
    os.writeEntryIfDifferent<scalar>("relax", 0, relax_);

    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<scalar>("M", 0, M_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletMachNumberPressureFvPatchScalarField
    );
}

// ************************************************************************* //
