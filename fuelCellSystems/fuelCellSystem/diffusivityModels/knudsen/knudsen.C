/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "knudsen.H"
#include "addToRunTimeSelectionTable.H"

#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusivityModels
{
    defineTypeNameAndDebug(knudsen, 0);
    addToRunTimeSelectionTable(diffusivityModel, knudsen, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

knudsen::knudsen
(
    const fvMesh& mesh,
    scalarField& diff,
    const labelList& cells,
    const dictionary& dict
)
:
    diffusivityModel(mesh, diff, cells, dict),
    Tname_(dict_.lookup("Tname")),
    dPore_(dict_.lookup("dPore")),
    MW_(dict_.lookup("MW"))
{}


knudsen::knudsen
(
    const fvMesh& mesh,
    scalarField& diff,
    const labelList& cells,
    word Tname,
    const dimensionedScalar& dPore,
    const dimensionedScalar& MW
)
:
    diffusivityModel(mesh, diff, cells),
    Tname_(Tname),
    dPore_(dPore),
    MW_(MW)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void knudsen::writeData()
{
    Info<< "diffusivityModels::knudsen:" << nl
        << "dPore = " << dPore_ << nl
        << "MW = " << MW_ << endl;
}


void knudsen::evaluate()
{
    //  D_{knudsen} = (poreDiameter/2)*97*sqrt(T/MW)
    //  where
    //      poreDiameter = [m]
    //      T ............ [K]
    //      MW ........... [kg/kmol]
    //  Geankoplis, Christie J, Transport Processes and Unit Operations,
    //  second edition (1983), Allyn and Bacon Series in Engineering,
    //  ISBN 0-205-07788-9, page 452.

#ifdef OF_VER_15
    // When using OpenFOAM-1.5, 1.5-dev
    const volScalarField& T =
        mesh_.db().lookupObject<volScalarField>(Tname_);
#else
    // When using OpenFOAM-1.6.x, ...
    const volScalarField& T =
        mesh_.thisDb().lookupObject<volScalarField>(Tname_);
#endif

    forAll(cells_, i)
    {
        diff_[cells_[i]] =
            48.5*dPore_.value()
            *Foam::sqrt(T.internalField()[cells_[i]]/MW_.value());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusivityModels
} // End namespace Foam

// ************************************************************************* //

