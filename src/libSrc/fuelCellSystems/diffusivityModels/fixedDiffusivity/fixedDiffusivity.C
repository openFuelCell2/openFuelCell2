/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "fixedDiffusivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace diffusivityModels
    {
        defineTypeNameAndDebug(fixedDiffusivity, 0);
        addToRunTimeSelectionTable(diffusivityModel, fixedDiffusivity, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModels::fixedDiffusivity::fixedDiffusivity
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
    const dictionary& dict
)
:
    diffusivityModel(name, mesh, diff, dict),
    diff0_(dict_.get<dimensionedScalar>("diff0"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::diffusivityModels::fixedDiffusivity::writeData()
{
    if (firstIndex_)
    {
        Info<< "diffusivityModels::fixedDiffusivity: " << nl
            << "   diff0 = " << diff0_ << endl;

        firstIndex_ = false;
    }
}


void Foam::diffusivityModels::fixedDiffusivity::evaluate()
{
    forAll(cells_, i)
    {
        diff_[cells_[i]] = diff0_.value();
    }
}

// ************************************************************************* //

