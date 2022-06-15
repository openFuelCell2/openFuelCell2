/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "diffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<diffusivityModel> diffusivityModel::New
(
    const word& name,
    const fvMesh& mesh,
    scalarField& diff,
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

    return autoPtr<diffusivityModel>
    (
        cstrIter()
        (
            name,
            mesh,
            diff,
            dict.subDict(diffTypeName + "Coeffs")
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
