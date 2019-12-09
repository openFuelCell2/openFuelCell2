/*---------------------------------------------------------------------------*\
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

    // Standard atmosphere pressure [Pa]
    // BIPM 10th Conferance Generale des Poids et Mesures (resolution 4)
    const dimensionedScalar pAtm("pAtm", dimPressure, 1.01325e5);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

binaryFSG::binaryFSG
(
    const fvMesh& mesh,
    scalarField& diff,
    const labelList& cells,
    const dictionary& dict
)
:
    diffusivityModel(mesh, diff, cells, dict),
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


binaryFSG::binaryFSG
(
    const fvMesh& mesh,
    scalarField& diff,
    const labelList& cells,
    word Tname,
    word pName,
    word spA,
    word spB
)
:
    diffusivityModel(mesh, diff, cells),
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


void binaryFSG::setSpecies(word spA, word spB)
{
    spA_ = spA;
    spB_ = spB;

    mA_ = fsgMolecularWeights(spA_);
    vA_ = fsgDiffusionVolumes(spA_);
 
    mB_ = fsgMolecularWeights(spB_);
    vB_ = fsgDiffusionVolumes(spB_);
}


void binaryFSG::writeData()
{
    Info<< "diffusivityModels::binaryFSG:" << nl;
    Info<< "    diffusing speciesA molWt diffVol: "
        << spA_ << " " << mA_ << " " << vA_ << nl;
    Info<< "   background speciesB molWt diffVol: "
        << spB_ << " " << mB_ << " " << vB_ << nl;
}


void binaryFSG::evaluate()
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

#ifdef OF_VER_15
    // When using OpenFOAM-1.5, 1.5-dev
    const volScalarField& T = mesh_.db().lookupObject<volScalarField>(Tname_);
    const volScalarField& p = mesh_.db().lookupObject<volScalarField>(pName_);
#else
    // When using OpenFOAM-1.6.x, ...
    const volScalarField& T =
        mesh_.thisDb().lookupObject<volScalarField>(Tname_);
    const volScalarField& p =
        mesh_.thisDb().lookupObject<volScalarField>(pName_);
#endif


    scalarField pTot = p/pAtm;

    forAll(cells_, i)
    {
        diff_[cells_[i]] =
            1e-7*pow(T[cells_[i]], 1.75)*sqrt(1/mA_ + 1/mB_)
            /(pTot[cells_[i]]*sqr(cbrt(vA_) + cbrt(vB_)));

    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// ************************************************************************* //

