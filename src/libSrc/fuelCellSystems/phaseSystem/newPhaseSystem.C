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

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseSystem> Foam::phaseSystem::New
(
    const fvMesh& mesh
)
{
    const word phaseSystemType
    (
        IOdictionary
        (
            IOobject
            (
                propertiesName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("type")
    );

    Info << "Selecting phaseSystemType for region "
        << mesh.name()
        << ": "
        << phaseSystemType << endl;

    auto* ctorPtr =
        dictionaryConstructorTable(phaseSystemType);

    if (!ctorPtr)
    {
        FatalErrorIn("phaseSystemType::New")
           << "Unknown phaseSystemType type "
           << phaseSystemType << endl << endl
           << "Valid phaseSystemType are : " << endl
           << dictionaryConstructorTablePtr_->sortedToc()
           << exit(FatalError);
    }

    return ctorPtr(mesh);
}


// ************************************************************************* //
