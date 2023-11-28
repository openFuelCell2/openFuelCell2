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

#include "porousFSG.H"
#include "addToRunTimeSelectionTable.H"

#include "labelList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace diffusivityModels
    {
        defineTypeNameAndDebug(porousFSG, 0);
        addToRunTimeSelectionTable(diffusivityModel, porousFSG, dictionary);
    }
}

// Standard atmosphere pressure [Pa]
// BIPM 10th Conferance Generale des Poids et Mesures (resolution 4)
const Foam::dimensionedScalar pAtm("pAtm", Foam::dimPressure, 1.01325e5);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::porousFSG::porousFSG
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
    alphaName_(dict_.lookup("alphaName")),
    spA_(dict_.lookup("speciesA")),
    spB_(dict_.lookup("speciesB")),
    eps_(1),
    tau_(1),
    dPore_(dict_.lookup("dPore"))
{
    // molecular weights and diffusion volumes
    // spA
    if(fsgMolecularWeights.found(spA_))
    {
        mA_ = fsgMolecularWeights(spA_);
        vA_ = fsgDiffusionVolumes(spA_);
    }
    else
    {
        FatalErrorIn
        (   "porousFSG::porousFSG(fvMesh& scalarField&, "
            "const labelList&, const dictionary&)"
        )   << "invalid species for FSG: " << spA_ << nl
            << "Valid species are :" << endl
            << fsgMolecularWeights.toc()
            << exit(FatalError);
    }
 
    // spB
    if(fsgMolecularWeights.found(spB_))
    {
        mB_ = fsgMolecularWeights(spB_);
        vB_ = fsgDiffusionVolumes(spB_);
    }
    else
    {
        FatalErrorIn
        (   "porousFSG::porousFSG(fvMesh& scalarField&, "
            "const labelList&, const dictionary&)"
        )   << "invalid species for FSG: " << spB_ << nl
            << "Valid species are :" << endl
            << fsgMolecularWeights.toc()
            << exit(FatalError);
    }

    // porosity
    if(dict_.readIfPresent("porosity", eps_))
    {
        if (eps_ <= 0.0 || eps_ > 1.0)
        {
            FatalIOErrorIn
            (
                "Foam::porousFSG::porousFSG"
                "(const fvMesh&, const scalarField&, const labelList&, "
                "const dictionary&)",
                dict_
            )
                << "out-of-range porosity value " << eps_
                << exit(FatalIOError);
        }
    }

    // tortuosity
    if(dict_.readIfPresent("tortuosity", tau_))
    {
        if (tau_ < 1.0)
        {
            FatalIOErrorIn
            (
                "Foam::porousFSG::porousFSG"
                "(const fvMesh&, const scalarField&, const labelList&, "
                "const dictionary&)",
                dict_
            )
                << "out-of-range tortuosity value " << tau_
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusivityModels::porousFSG::setSpecies(word spA, word spB)
{
    spA_ = spA;
    spB_ = spB;

    mA_ = fsgMolecularWeights(spA_);
    vA_ = fsgDiffusionVolumes(spA_);

    mB_ = fsgMolecularWeights(spB_);
    vB_ = fsgDiffusionVolumes(spB_);
}



void Foam::diffusivityModels::porousFSG::writeData()
{
    if (firstIndex_)
    {
        Info<< "diffusivityModels::porousFSG:" << nl;
        Info<< "    diffusing speciesA molWt diffVol:  "
            << spA_ << " " << mA_ << " " << vA_ << nl;
        Info<< "    background speciesB molWt diffVol: "
            << spB_ << " " << mB_ << " " << vB_ << nl;
        Info<< "    dPore      = " << dPore_ << nl
            << "    porosity   = " << eps_ << nl
            << "    tortuosity = " << tau_ << endl;

        firstIndex_ = false;
    }
}


void Foam::diffusivityModels::porousFSG::evaluate()
{
    // binaryFSG
    // ---------
    //if (doBinary_)
    //{
        binaryFSG bfsg =
               binaryFSG(this->zoneName_, this->mesh_, this->diff_,
                         this->Tname_, this->pName_,
                         this->spA_, this->spB_);
        bfsg.evaluate();
    //}

    // knudsen
    // -------
    scalarField diffK(diff_);
    dimensionedScalar MW ("MW", dimensionSet(1,0,0,0,-1,0,0), mA_);

    knudsen knud =
        knudsen(this->zoneName_, this->mesh_, diffK,
                this->Tname_, this->dPore_, MW);
    knud.evaluate();

    //- alpha value
    const volScalarField& alpha =
        this->mesh_.lookupObject<volScalarField>(alphaName_); 

    // combine as (eps/tau)*harmonicMean
    // -------
    forAll(cells_, i)
    {
        diff_[cells_[i]] =
            Foam::pow(alpha[cells_[i]], 1.5)
          * eps_/tau_/(1/diff_[cells_[i]] + 1/diffK[cells_[i]]);
            //eps_*diff_[cells_[i]]*diffK[cells_[i]]/
            //(tau_*(diff_[cells_[i]] + diffK[cells_[i]]));
    }
}

// ************************************************************************* //

