/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    
\*---------------------------------------------------------------------------*/

#include "binaryFSG.H"
#include "addToRunTimeSelectionTable.H"

#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace diffusivityModels
    {
        defineTypeNameAndDebug(binaryFSG, 0);
        addToRunTimeSelectionTable(diffusivityModel, binaryFSG, dictionary);
    }
}

// Standard atmosphere pressure [Pa]
// BIPM 10th Conferance Generale des Poids et Mesures (resolution 4)
const Foam::dimensionedScalar pAtm("pAtm", Foam::dimPressure, 1.01325e5);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::binaryFSG::binaryFSG
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
    const dictionary& dict
)
:
    diffusivityModel(name, mesh, diff, dict),
    Tname_(dict_.lookup("Tname")),
    pName_(dict_.lookup("pName")),
    spA_(dict_.lookup("speciesA")),
    spB_(dict_.lookup("speciesB"))
{
    if(fsgMolecularWeights.found(spA_))
    {
        mA_ = fsgMolecularWeights(spA_);
        vA_ = fsgDiffusionVolumes(spA_);
    }
    else
    { 
        FatalErrorIn
        (   "binaryFSG::binaryFSG(fvMesh& scalarField&, const labelList&, "
            "const dictionary&)"
        )   << "invalid species for FSG: " << spA_ << nl
            << "Valid species are :" << endl
            << fsgMolecularWeights.toc()
            << exit(FatalError);
    }
 
    if(fsgMolecularWeights.found(spB_))
    {
        mB_ = fsgMolecularWeights(spB_);
        vB_ = fsgDiffusionVolumes(spB_);
    }
    else
    { 
        FatalErrorIn
        (   "binaryFSG::binaryFSG(fvMesh& scalarField&, const labelList&, "
            "const dictionary&)"
        )   << "invalid species for FSG: " << spB_ << nl
            << "Valid species are :" << endl
            << fsgMolecularWeights.toc()
            << exit(FatalError);
    }
}


Foam::diffusivityModels::binaryFSG::binaryFSG
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
    word Tname,
    word pName,
    word spA,
    word spB
)
:
    diffusivityModel(name, mesh, diff),
    Tname_(Tname),
    pName_(pName),
    spA_(spA),
    spB_(spB)
{
    mA_ = fsgMolecularWeights(spA_);
    vA_ = fsgDiffusionVolumes(spA_);
 
    mB_ = fsgMolecularWeights(spB_);
    vB_ = fsgDiffusionVolumes(spB_);
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusivityModels::binaryFSG::setSpecies(word spA, word spB)
{
    spA_ = spA;
    spB_ = spB;

    mA_ = fsgMolecularWeights(spA_);
    vA_ = fsgDiffusionVolumes(spA_);

    mB_ = fsgMolecularWeights(spB_);
    vB_ = fsgDiffusionVolumes(spB_);
}


void Foam::diffusivityModels::binaryFSG::writeData()
{
    if (firstIndex_)
    {
        Info<< "diffusivityModels::binaryFSG:" << nl;
        Info<< "    diffusing speciesA molWt diffVol: "
            << spA_ << " " << mA_ << " " << vA_ << nl;
        Info<< "   background speciesB molWt diffVol: "
            << spB_ << " " << mB_ << " " << vB_ << nl;

        firstIndex_ = false;
    }
}


void Foam::diffusivityModels::binaryFSG::evaluate()
{
    //             1e-3 * T^{1.75} * sqrt(1/mA + 1/mB)
    //  D = 1e-4 * -----------------------------------
    //                p * [ vA^(1/3) + vB^{1/3} ]^2
    //  where
    //      D = diffusivity ......... [m^2/s]
    //      T = temperature ......... [K]
    //      p = total pressure ...... [atm]
    //      m = molecular weight .... [kg/kmol]
    //      v = diffusion volume .... [cm^3]     NOTE: cm
    //      A,B = species index
    //
    //  Fuller, Schettler, and Giddings,
    //  A new method for prediction of binary gas-phase diffusion coefficients,
    //  Industrial and Engineering Chemistry, v58, n5, May, 1966, pp 19-27.

    using Foam::pow;
    using Foam::sqr;
    using Foam::sqrt;
    using Foam::cbrt;

    const volScalarField& T =
        mesh_.lookupObject<volScalarField>(Tname_);
    const volScalarField& p =
        mesh_.lookupObject<volScalarField>(pName_);


    scalarField pTot = p/pAtm;

    forAll(cells_, i)
    {
        diff_[cells_[i]] =
            1e-7*pow(T[cells_[i]], 1.75)*sqrt(1/mA_ + 1/mB_)
            /(pTot[cells_[i]]*sqr(cbrt(vA_) + cbrt(vB_)));

    }
}

// ************************************************************************* //

