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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDiffusivity::fixedDiffusivity
(
    const fvMesh& mesh,
    scalarField& diff,
    const labelList& cells,
    const dictionary& dict
)
:
    diffusivityModel(mesh, diff, cells, dict),
    diff0_(dict_.lookup("diff0"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void fixedDiffusivity::writeData()
{
    Info<< "diffusivityModels::fixedDiffusivity: " << nl
        << "   diff0 = " << diff0_ << endl;
}


void fixedDiffusivity::evaluate()
{
    forAll(cells_, i)
    {
	diff_[cells_[i]] = diff0_.value();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// ************************************************************************* //

