/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<diffusivityModel> diffusivityModel::New
(
    const fvMesh& mesh,
    scalarField& diff,
    const labelList& cells,
    const dictionary& dict
)
{
    word diffTypeName = dict.get<word>("type");

    Info<< "Selecting diffusivity model " << diffTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(diffTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "diffusivityModel::New(fvMesh& mesh, scalarField&, "
            "const labelList&)"
        )   << "Unknown diffusivityModel type " << diffTypeName
            << endl << endl
            << "Valid diffusivityModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<diffusivityModel>(cstrIter()(mesh, diff, cells, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// ************************************************************************* //
