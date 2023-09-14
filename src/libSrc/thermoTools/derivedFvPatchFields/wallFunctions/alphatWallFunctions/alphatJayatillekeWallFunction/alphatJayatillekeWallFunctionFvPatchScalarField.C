/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd
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

#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphatJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


tmp<scalarField> alphatJayatillekeWallFunctionFvPatchScalarField::yPlus
(
    const compressible::turbulenceModel& turbModel
) const
{
    const label patchi = patch().index();

    const tmp<volScalarField> tnut = turbModel.nut();
    const volScalarField& nut = tnut();

    if (isA<nutWallFunctionFvPatchScalarField>(nut.boundaryField()[patchi]))
    {
        const auto& nutPf =
            dynamic_cast<const nutWallFunctionFvPatchScalarField&>
            (
                nut.boundaryField()[patchi]
            );

        return nutPf.yPlus();
    }
    else
    {
        const scalarField& y = turbModel.y()[patchi];
        const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

        return
            y*sqrt(turbModel.nuEff(patchi)*mag(Uw.snGrad()))
           /turbModel.nu(patchi);
    }
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11;

    for (int iter = 0; iter < maxIters_; ++iter)
    {
        const scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        const scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        const scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.getOrDefault<scalar>("Prt", 0.85)),
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    E_(dict.getOrDefault<scalar>("E", 9.8))
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField yPlusp(yPlus(turbModel));

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    const scalarField& alphaw = talphaw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& hew =
        turbModel.transport().he().boundaryField()[patchi];

    scalarField& alphatw = *this;

    // Heat flux [W/m2] - lagging alphatw
    const scalarField qDot
    (
        turbModel.transport().alphaEff(alphatw, patchi)*hew.snGrad()
    );

    // Populate boundary values
    forAll(alphatw, facei)
    {
        const scalar yPlus = yPlusp[facei];

        const scalar uTau = yPlus/y[facei]*(muw[facei]/rhow[facei]);

        // Molecular Prandtl number
        const scalar Pr = muw[facei]/alphaw[facei];

        // Molecular-to-turbulent Prandtl number ratio
        const scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        const scalar P = Psmooth(Prat);
        const scalar yPlusTherm = this->yPlusTherm(P, Prat);

        // Evaluate new effective thermal diffusivity
        scalar alphaEff = 0;
        if (yPlus < yPlusTherm)
        {
            const scalar A = qDot[facei]*rhow[facei]*uTau*y[facei];
            const scalar B = qDot[facei]*Pr*yPlus;
            const scalar C = Pr*0.5*rhow[facei]*uTau*sqr(magUp[facei]);

            alphaEff = A/(B + C + VSMALL);
        }
        else
        {
            const scalar A = qDot[facei]*rhow[facei]*uTau*y[facei];
            const scalar B = qDot[facei]*Prt_*(1.0/kappa_*log(E_*yPlus) + P);
            const scalar magUc =
                uTau/kappa_*log(E_*yPlusTherm) - mag(Uw[facei]);
            const scalar C =
                0.5*rhow[facei]*uTau
               *(Prt_*sqr(magUp[facei]) + (Pr - Prt_)*sqr(magUc));

            alphaEff = A/(B + C + VSMALL);
        }

        // Update turbulent thermal diffusivity
        alphatw[facei] = max(scalar(0), alphaEff - alphaw[facei]);

        if (debug)
        {
            Info<< "    uTau           = " << uTau << nl
                << "    Pr             = " << Pr << nl
                << "    Prt            = " << Prt_ << nl
                << "    qDot           = " << qDot[facei] << nl
                << "    yPlus          = " << yPlus << nl
                << "    yPlusTherm     = " << yPlusTherm << nl
                << "    alphaEff       = " << alphaEff << nl
                << "    alphaw         = " << alphaw[facei] << nl
                << "    alphatw        = " << alphatw[facei] << nl
                << endl;
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<scalar>("Prt", 0.85, Prt_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("E", 9.8, E_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
