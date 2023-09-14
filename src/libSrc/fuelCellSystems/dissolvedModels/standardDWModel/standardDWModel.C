/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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

#include "fvCFD.H"
#include "standardDWModel.H"
#include "regionType.H"
#include "fuelCellSystem.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dissolvedModels
{
    defineTypeNameAndDebug(standardDWModel, 0);

    addToRunTimeSelectionTable
    (
        dissolvedModel,
        standardDWModel,
        dictionary
    );
}
}

const Foam::dimensionedScalar F = Foam::constant::physicoChemical::F;

// * * * * * * * * * * * * * * * * Private functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::dissolvedModels::standardDWModel::DNaf
(
    const scalar& d0,
    const scalar& w
) const
{
    scalar tmp(0.0);

    if(w <= 3.0)
    {
        tmp = d0*3.1e-7*w*(pow(M_E, 0.28*w) - 1.0);
    }
    else
    {
        tmp = d0*4.17e-8*w*(161*pow(M_E, -w) + 1.0);
    }

    return tmp;
}


Foam::scalar Foam::dissolvedModels::standardDWModel::effLam(const scalar& act) const
{
    scalar tmp(0.0);

    if(act <= 1.0)
    {
        tmp = 0.043 + 17.81*act - 39.85*Foam::pow(act, 2) + 36.0*Foam::pow(act, 3);
    }
    else if (act > 3.0)
    {
        tmp = 15.8;
    }
    else
    {
        tmp = 14.0 + 1.4*(act - 1.0);
    }

    return tmp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissolvedModels::standardDWModel::standardDWModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dissolvedModel(mesh),

    dict_(dict),
    ksi_(dict_.lookup("ksi")),
    nd_("nd", dimless, dict_),
    rhoOnEW_("rhoOnEW", dimMoles/dimVol, dict_),
    iName_(dict_.lookupOrDefault<word>("i", "i")),
    TName_(dict_.lookupOrDefault<word>("T", "T")),
    relax_(dict_.lookupOrDefault<scalar>("relax", 0.0)),
    corr_(dict_.lookupOrDefault<scalar>("corr", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dissolvedModels::standardDWModel::~standardDWModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::dissolvedModels::standardDWModel::solve()
{
    const volVectorField& i = mesh().lookupObject<volVectorField>(iName_);

    surfaceScalarField phiI = fvc::interpolate(i) & mesh().Sf();

    tmp<fvScalarMatrix> wEqn
    (
        fvm::ddt(rhoOnEW_, lambda_)
      + fvm::div(phiI*nd_/F, lambda_, "div(D,lambda)")
      - fvm::laplacian(rhoOnEW_*Dwm_, lambda_, "laplacian(diff,lambda)")
      - dmdt_
    );

    //- Set reference values
    {
        const scalarField& source = dmdt_;
        const scalarField& volume = mesh().V();

        // Update the water content
        // Water intake = water uptake
        scalarField sum = source*volume;
        scalar iDot(Foam::gSum(sum)/Foam::gSum(volume));

        for (label i = 0; i < lambda_.size(); i++)
        {
            wEqn->setReference(i, lambda_[i] + iDot*relax_);
        }
    }

    wEqn->relax();
    wEqn->solve();
}


void Foam::dissolvedModels::standardDWModel::correct()
{
    scalarField& lambdaIn = lambda_;

    const scalarField& T = mesh().lookupObject<volScalarField>(TName_);

    scalarField Dwm0 = Foam::exp(-2436/T)*corr_;

    forAll(Dwm_, cellI)
    {
        Dwm_[cellI] = DNaf(Dwm0[cellI], lambdaIn[cellI]);
    }

    Dwm_.correctBoundaryConditions();
    act_.correctBoundaryConditions();
}


void Foam::dissolvedModels::standardDWModel::update
(
    const word& clName
)
{
    scalarField& lambdaIn = lambda_;

    //- only consider the catalyst layer for water adsorption or desorption

    label znId = mesh().cellZones().findZoneID(clName);
    const labelList& cells = mesh().cellZones()[znId];

    forAll(cells, cellI)
    {
        label cellId = cells[cellI];

        //- dmdt = ksi*rhoOnEW*(lambda^eff - lambda)
        dmdt_[cellId] =
            ksi_[clName]*rhoOnEW_.value()
          * (effLam(act_[cellId]) - lambdaIn[cellId]);
    }

    dmdt_.correctBoundaryConditions();
}


void Foam::dissolvedModels::standardDWModel::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    // nothing
}


void Foam::dissolvedModels::standardDWModel::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    // nothing
}


bool Foam::dissolvedModels::standardDWModel::read(const dictionary& dict)
{
    const dictionary& dict0 = dict.subDict(type() + "Coeffs");

    dict0.lookup("ksi") >> ksi_;
    nd_ = dimensionedScalar("nd", dimless, dict0);
    rhoOnEW_ = dimensionedScalar("rhoOnEW", dimMoles/dimVol, dict0);
    iName_ = word(dict0.lookupOrDefault<word>("i", "i"));
    TName_ = word(dict0.lookupOrDefault<word>("T", "T"));
    relax_ = scalar(dict0.lookupOrDefault<scalar>("relax", 0.0));
    corr_ = scalar(dict0.lookupOrDefault<scalar>("corr", 1.0));

    return true;
}
// ************************************************************************* //
